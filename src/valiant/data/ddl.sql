/*
EXPERIMENT-wide
*/

-- Targetons

create table targetons (
    id integer primary key,
    start integer not null,
    end integer not null,
    chr text not null,
    strand text not null,

    check (
        start >= 1 and
        end >= start and
        (strand in ('+', '-'))
    )
);

-- Transcript exons

create table exons (
    id integer primary key,
    start integer not null,
    end integer not null,
    exon_index integer not null
);
create unique index exons_exon_index_idx on exons (exon_index);

-- Compute the CDS prefix lengths
create view if not exists v_exon_ext as
select
    id,
    start,
    end, (
        lead((end - start + 1) % 3, -1, 0)
        over (order by exon_index)
    ) as cds_prefix_length
from exons;

-- Background edits

create table background_variants (
    id integer primary key,
    var_id text,
    start integer not null,
    -- start_bg integer not null,
    ref text,
    alt text
);

create view v_background_variants as
select *, alt_len - ref_len as alt_ref_delta
from (
    select
        start,
        ref,
        alt,
        length(ref) as ref_len,
        length(alt) as alt_len
    from background_variants
);

create table targeton_background_variants (
    targeton_id integer not null primary key references targetons (id),
    bg_id integer not null references background_variants (id),
    start integer not null  -- cell-line position
);

-- PAM protection edits

create table pam_protection_edits (
    id integer primary key,
    start integer not null,
    ref text not null,
    alt text not null
);

create table sgrna_ids (
    id integer primary key,
    name text not null
);
create unique index sgrna_ids_name_idx on sgrna_ids (name);

create table pam_protection_edit_sgrna_ids (
    var_ppe_id integer not null primary key references pam_protection_edits (id),
    sgrna_id integer not null references sgrna_ids (id)
) without rowid;

create table exon_codon_ppes (
    exon_id integer not null references exons (id),
    codon_index integer not null,
    -- Assumption: at most one PPE per codon
    ppe_id integer not null references pam_protection_edits (id),

    primary key (exon_id, codon_index)
);

-- Custom variants

create table custom_variant_collections (
    id integer primary key,
    name text not null
);

create table custom_variants (
    id integer primary key,
    var_id text,
    collection_id integer not null references custom_variant_collections (id),
    start integer not null,
    vcf_nt text not null,
    ref text not null,
    alt text not null
);

/*
PER TARGETON
*/

-- Pattern variants

create table pattern_variants (
    id integer primary key,
    mutator text not null,

    -- Altered reference
    pos_a integer not null,
    ref_a text not null,
    alt_a text not null,

    -- Oligonucleotide sequence
    oligo text not null,

    -- CDS-specific annotation
    -- Possibly partial codon including the PPE, if any
    codon_ref_a text,
    codon_alt_a text,
    -- Amino acid change
    aa_ref text,
    aa_alt text
);

-- Map positions to reference
create view if not exists v_pattern_variants as
select
    *
from pattern_variants pv;

create table targeton_pam_protection_edits (
    id integer primary key references pam_protection_edits (id),
    start integer not null
);

create table targeton_custom_variants (
    id integer primary key references custom_variants (id),
    start integer not null,
    in_const int not null default 0
);

create table mutations (
    id integer primary key,
    -- Original reference
    pos_r integer not null,
    ref_r text not null,
    alt_r text not null,
    -- Altered reference
    pos_a integer not null,
    ref_a text not null,
    alt_a text not null,
    -- CDS-specific annotation
    codon_ref text,
    codon_alt text,
    aa_ref text,
    aa_alt text
);

create view if not exists v_meta_custom as
select
    cv.start as ref_start,
    cv.ref,
    cv.alt,
    null as ref_aa,
    null as alt_aa,
    cv.var_id as vcf_var_id,
    cvc.name as vcf_alias,
    'custom' as mutator,
    tcv.in_const
from targeton_custom_variants tcv
left join custom_variants cv on cv.id = tcv.id
left join custom_variant_collections cvc on cvc.id = cv.collection_id;

create view if not exists v_meta as
select
    s.ref_start,
    s.ref,
    s.alt,
    s.ref_aa,
    s.alt_aa,
    s.vcf_var_id,
    s.vcf_alias,
    s.mutator,
    s.in_const,
    -- TODO: handle multiple guides
    si.name as sgrna_ids
from (
    select
        w.*,
        e.id as exon_id,
        ((w.ref_start - e.start - e.cds_prefix_length) / 3) as codon_index
    from (
        select
            v.*,
            v.ref_start + min(0, length(v.ref) - 1) as ref_end
        from (
            select
                pos_a as ref_start,  -- pos_r?
                ref_a as ref,
                alt_a as alt,
                aa_ref as ref_aa,
                aa_alt as alt_aa,
                null as vcf_var_id,
                null as vcf_alias,
                mutator,
                0 as in_const
            from pattern_variants pv
            union all
            select * from v_meta_custom
        ) v
    ) w
    left join v_exon_ext e on (
        (w.ref_start >= e.start and w.ref_start <= e.end) or
        (w.ref_end >= e.start and w.ref_end <= e.end)
    )
) s
left join exon_codon_ppes ecp on
    ecp.exon_id = s.exon_id and
    ecp.codon_index = s.codon_index
left join pam_protection_edit_sgrna_ids ppesi on
    ppesi.var_ppe_id = ecp.ppe_id
left join sgrna_ids si on si.id = ppesi.sgrna_id;
