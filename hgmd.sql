create table annovcf (
    CHROM char(2),
    POS int,
    REF text,
    ALT text,
    QUAL numeric,
    FILTER text,
    AC int,
    AN int,
    AF numeric,
    AA_CHANGE text,
    HGMD_MUT text,
    HGMD_GENE text,
    HGMD_SITE text,
    HGNC_GENE text,
    LOF text,
    LOF_FLAG text
);

create table masites (
    CHROM char(2),
    POS int
);

# table allmut was created from a mysql dump: /humgen/atgu1/fs03/eminikel/057hgmd/allmut.sql

select   *
from     allmut
where    tag = 'DM'
;

select   *
from     allmut a, masites ma, annovcf v
where    a.chromosome = v.CHROM
and      a.startCoord = v.POS
and      