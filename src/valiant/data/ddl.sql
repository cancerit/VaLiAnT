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

-- Gene offsetting (just transient?)

create table background_offsets (
    ref_pos integer not null primary key,
    offset integer not null
);

-- Transcript exons

create table exons (
    id integer primary key,
    start integer not null,
    end integer not null,
    exon_index integer not null
);
create unique index exons_exon_index_idx on exons (exon_index);

-- Background edits

create table background_variants (
    id integer primary key,
    var_id text,
    start integer not null,
    -- start_bg integer not null,
    ref text,
    alt text
);

create view background_variants_v as
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

-- Custom variants

create table custom_variant_collections (
    id integer primary key,
    name text not null
);

-- Include generated variants? How to express the source in either (mutator vs. collection)?
create table custom_variants (
    id integer primary key,
    var_id text,
    collection_id integer not null references custom_variant_collections (id),
    start integer not null,
    vcf_nt text not null,
    ref text not null,
    alt text not null
);

-- Pattern variants

create table pattern_variants (
    id integer primary key,
    mutator text not null,
    -- Original reference
    pos_r integer,
    ref_r text,
    alt_r text,
    -- Altered reference
    pos_a integer not null,
    ref_a text not null,
    alt_a text not null,

    -- CDS-specific annotation
    -- Possibly partial codon including the PPE, if any
    codon_ref_a text,
    codon_alt_a text,
    -- Amino acid change
    aa_ref text,
    aa_alt text
);
