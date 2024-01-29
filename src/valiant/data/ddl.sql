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
    exon_index integer not null,
    cds_prefix_length integer not null,
    cds_suffix_length integer not null,
    first_codon_start integer not null
);
create unique index exons_exon_index_idx on exons (exon_index);

-- Compute the CDS prefix lengths
create view v_exon_ext as
select
    e.*,
    e.length / 3 as last_codon_index,
    e.cds_prefix_length as frame  -- to avoid breaking other queries
from (
    select *, (end - start + 1) as length
    from exons
) e;

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
select
    *,
    alt_len - ref_len as alt_ref_delta,
    start + max(0, ref_len - 1) as ref_end
from (
    select
        id,
        start,
        ref,
        alt,
        length(ref) as ref_len,
        length(alt) as alt_len
    from background_variants
);

create view v_background_variant_codons as
select
    t.*,
    start_delta / 3 as start_codon_index,
    start_delta % 3 as start_codon_offset,
    end_delta / 3 as end_codon_index,
    end_delta % 3 as end_codon_offset
from (
    select
        b.*,
        -- Start codon
        es.exon_index as start_exon_index,
        es.last_codon_index as start_last_codon_index,
        abs(b.start - es.first_codon_start) as start_delta,
        -- End codon
        ee.exon_index as end_exon_index,
        ee.last_codon_index as end_last_codon_index,
        abs(b.ref_end - ee.first_codon_start) as end_delta
    from v_background_variants b
    left join v_exon_ext es on
        b.start >= es.start and
        b.start <= es.end
    left join v_exon_ext ee on
        b.ref_end >= ee.start and
        b.ref_end <= ee.end
) t;

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

create view v_ppe_sgrna_ids as
select ppe.*, si.name as sgrna_id
from pam_protection_edits as ppe
left join pam_protection_edit_sgrna_ids ppesi on ppesi.var_ppe_id = ppe.id
left join sgrna_ids si on si.id = ppesi.sgrna_id;

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

create view v_custom_variants as
select *, start + max(0, length(ref) - 1) as end
from custom_variants;

/*
PER TARGETON
*/

-- PAM protection edits

create table targeton_pam_protection_edits (
    id integer primary key references pam_protection_edits (id),
    start integer not null
);

create table targeton_exon_codon_ppes (
    exon_id integer not null references exons (id),
    codon_index integer not null,
    -- Assumption: at most one PPE per codon
    ppe_id integer not null references pam_protection_edits (id),

    primary key (exon_id, codon_index)
);

create view v_exon_codon_ppes as
select
    ecp.exon_id,
    e.exon_index,
    ecp.codon_index,
    ecp.ppe_id,
    ppe.start as ppe_start,
    abs(ppe.start - e.first_codon_start) % 3 as codon_offset
from targeton_exon_codon_ppes ecp
left join v_exon_ext e on e.id = ecp.exon_id
left join pam_protection_edits ppe on ppe.id = ecp.ppe_id;

create view v_exon_ppes as
select
    e.exon_index,
    ppe.start,
    ppe.ref,
    ppe.alt,
    si.name as sgrna_id
from targeton_pam_protection_edits tppe
left join pam_protection_edits ppe on ppe.id = tppe.id
left join pam_protection_edit_sgrna_ids ppesi on ppesi.var_ppe_id = ppe.id
left join sgrna_ids si on si.id = ppesi.sgrna_id
inner join exons e on
    ppe.start >= e.start and
    ppe.start <= e.end;

-- Pattern variants

create table alt_pattern_variants (
    id integer primary key,
    mutator text not null,

    -- Original reference
    pos_r integer not null,

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
    aa_alt text,
    mutation_type text  -- syn|mis|non
);

create table ref_pattern_variants (
    id integer primary key references alt_pattern_variants (id),
    start integer not null,
    ref text not null,
    alt text not null
);

create table targeton_custom_variants (
    id integer primary key references custom_variants (id),
    start integer not null,
    oligo text not null,
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

create view v_meta_pattern as
select
    pos_r as original_start,  -- to be computed
    pos_a as ref_start,
    ref_a as ref,
    alt_a as alt,
    aa_ref as ref_aa,
    aa_alt as alt_aa,
    null as vcf_var_id,
    null as vcf_alias,
    '' as vcf_nt,  -- to be computed
    mutator,
    0 as in_const,
    oligo,
    mutation_type
from alt_pattern_variants;

create view v_meta_custom as
select
    cv.start as original_start,
    tcv.start as ref_start,
    cv.ref,
    cv.alt,
    null as ref_aa,
    null as alt_aa,
    cv.var_id as vcf_var_id,
    cvc.name as vcf_alias,
    vcf_nt,
    'custom' as mutator,
    tcv.in_const,
    tcv.oligo,
    null as mutation_type
from targeton_custom_variants tcv
left join custom_variants cv on cv.id = tcv.id
left join custom_variant_collections cvc on cvc.id = cv.collection_id;

create view v_meta as
-- Map: PPE position -> sgRNA ID
with t_pos_sgrna_ids as (
    select
        ppe.start,
        si.name as sgrna_id
    from targeton_pam_protection_edits tppe
    left join pam_protection_edits ppe on ppe.id = tppe.id
    left join pam_protection_edit_sgrna_ids ppesi on
        ppesi.var_ppe_id = ppe.id
    left join sgrna_ids si on
        si.id = ppesi.sgrna_id
)
select
    s.original_start,
    s.ref_start,
    s.ref_end,
    s.ref,
    s.alt,
    s.ref_aa,
    s.alt_aa,
    s.vcf_var_id,
    s.vcf_alias,
    s.vcf_nt,
    s.mutator,
    s.in_const,
    s.oligo,
    s.mutation_type,
    s.start_exon_index,
    ecps.ppe_start as start_ppe_start,
    s.start_codon_index,
    s.end_exon_index,
    s.end_codon_index,
    ecpe.ppe_start as end_ppe_start,
    (
        select
            group_concat(z.sgrna_id, ';')
        from (
            -- Filter PPE's by position
            select sgrna_id
            from t_pos_sgrna_ids as t
            where
                t.start >= ifnull(ecps.ppe_start, s.ref_start) and
                t.start <= ifnull(ecpe.ppe_start, s.ref_end)
            group by t.sgrna_id
            order by t.sgrna_id
        ) z
    ) as sgrna_ids
from (
    select
        w.*,
        es.id as start_exon_id,
        es.exon_index as start_exon_index,
        es.last_codon_index as start_last_codon_index,
        (abs(w.ref_start - es.first_codon_start) / 3) as start_codon_index,
        ee.id as end_exon_id,
        ee.exon_index as end_exon_index,
        ee.last_codon_index as end_last_codon_index,
        (abs(w.ref_end - ee.first_codon_start) / 3) as end_codon_index
    from (
        select
            v.*,
            v.ref_start + max(0, length(v.ref) - 1) as ref_end
        from (
            select * from v_meta_pattern
            union all
            select * from v_meta_custom
        ) v
    ) w
    left join v_exon_ext es on
        w.ref_start >= es.start and
        w.ref_start <= es.end
    left join v_exon_ext ee on
        w.ref_end >= ee.start and
        w.ref_end <= ee.end
) s
left join v_exon_codon_ppes ecps on
    ecps.exon_index = s.start_exon_index and
    ecps.codon_index = s.start_codon_index
left join v_exon_codon_ppes ecpe on
    ecpe.exon_index = s.end_exon_index and
    ecpe.codon_index = s.end_codon_index;
