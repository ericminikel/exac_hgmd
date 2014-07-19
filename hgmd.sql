create table annovcf (
    CHROM char(2),
    POS int,
    REF text,
    ALT text,
    QUAL numeric,
    FILTER text,
    AC int,
    AN int,
    AF float,
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

create index v_chrom on annovcf(CHROM);
create index v_pos   on annovcf(POS);

create index a_chrom on allmut(chromosome);
create index a_pos   on allmut(startCoord);

create index m_chrom on masites(CHROM);
create index m_pos   on masites(POS);

select   count(*)
from     allmut a
;
-- 144576

-- how many mutations _of interest_?
select   count(*)
from     allmut
where    tag = 'DM'
and      not disease like '%?'
;
-- 127625

select   count(*)
from     annovcf v
;
-- 24224

select   count(*)
from     allmut a, annovcf v
where    a.chromosome = v.CHROM
and      a.startCoord = v.POS
;
-- 24591
-- some many:many matches here.

select   count(*)
from     masites m, annovcf v
where    m.CHROM = v.CHROM
and      m.POS = v.POS
;
-- 6947

select   count(*)
from     annovcf v
where    not exists (
             select null from masites ma where ma.CHROM = v.CHROM and ma.POS = v.POS
         )
;
-- 17277

select   count(*)
from     annovcf v, allmut a
where    a.chromosome = v.CHROM
and      a.startCoord = v.POS
and      not exists (
             select null from masites ma where ma.CHROM = v.CHROM and ma.POS = v.POS
         )
;
-- 17533
-- so even removing the ma sites, there is some multiplicity

-- what about joining on HGMD_MUT id?
select   count(*)
from     allmut a, annovcf v
where    a.acc_num = v.HGMD_MUT
;
-- 24112

-- is acc_num unique in allmut?
select   count(*)
from     allmut a
group by acc_num
having   count(*) > 1
;
-- empty set. yes.

-- is HGMD_MUT unique in annovcf?
select   count(*)
from     annovcf v
group by HGMD_MUT
having   count(*) > 1
order by count(*) desc
limit 10
;
-- no

-- how about after removing ma sites?
select   count(*)
from     annovcf v
where     not exists (
             select null from masites ma where ma.CHROM = v.CHROM and ma.POS = v.POS
         )
group by HGMD_MUT
having   count(*) > 1
order by count(*) desc
limit 10
;
-- empty set. yes, unique then.



select   count(*)
from     allmut a, annovcf v
where    a.acc_num = v.HGMD_MUT
and      a.tag = 'DM'
and      not a.disease like '%?%'
and      not exists (
             select null from masites ma where ma.CHROM = v.CHROM and ma.POS = v.POS
         )
;
-- 17149

select   a.disease, a.gene, v.AF
from     allmut a, annovcf v
where    a.acc_num = v.HGMD_MUT
and      a.tag = 'DM'
and      not a.disease like '%?%'
and      not exists (
             select null from masites ma where ma.CHROM = v.CHROM and ma.POS = v.POS
         )
order by v.AF desc
limit 10
;

select   a.gene, v.AF, v.CHROM, v.POS, v.REF, v.ALT, a.pmid
from     allmut a, annovcf v
where    a.acc_num = v.HGMD_MUT
and      a.tag = 'DM'
and      not a.disease like '%?%'
and      not exists (
             select null from masites ma where ma.CHROM = v.CHROM and ma.POS = v.POS
         )
order by v.AF desc
limit 10
;
