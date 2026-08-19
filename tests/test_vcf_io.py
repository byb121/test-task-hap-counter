import textwrap

from hap_counter.vcf_io import SnvSite, read_biallelic_snvs

VCF_HEADER = textwrap.dedent(
    """\
    ##fileformat=VCFv4.2
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##FILTER=<ID=LowQual,Description="Low quality variant">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##contig=<ID=chr1,length=1000000>
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample
    """
)


def write_vcf(tmp_path, records):
    path = tmp_path / "variants.vcf"
    path.write_text(VCF_HEADER + "\n".join(records) + ("\n" if records else ""))
    return str(path)


def test_keeps_pass_biallelic_snv(tmp_path):
    vcf_path = write_vcf(
        tmp_path,
        ["chr1\t100\t.\tG\tA\t30\tPASS\t.\tGT\t1|0"],
    )
    assert list(read_biallelic_snvs(vcf_path)) == [SnvSite("chr1", 100, "G", "A")]


def test_skips_non_pass_filter(tmp_path):
    vcf_path = write_vcf(
        tmp_path,
        ["chr1\t100\t.\tG\tA\t30\tLowQual\t.\tGT\t1|0"],
    )
    assert list(read_biallelic_snvs(vcf_path)) == []


def test_skips_indel(tmp_path):
    vcf_path = write_vcf(
        tmp_path,
        ["chr1\t100\t.\tCAGAAT\tC\t30\tPASS\t.\tGT\t0/1"],
    )
    assert list(read_biallelic_snvs(vcf_path)) == []


def test_skips_multiallelic(tmp_path):
    vcf_path = write_vcf(
        tmp_path,
        ["chr1\t100\t.\tG\tA,T\t30\tPASS\t.\tGT\t1/2"],
    )
    assert list(read_biallelic_snvs(vcf_path)) == []


def test_skips_symbolic_allele(tmp_path):
    vcf_path = write_vcf(
        tmp_path,
        ["chr1\t100\t.\tG\t*\t30\tPASS\t.\tGT\t0/1"],
    )
    assert list(read_biallelic_snvs(vcf_path)) == []


def test_keeps_unphased_and_hom_alt_sites(tmp_path):
    vcf_path = write_vcf(
        tmp_path,
        [
            "chr1\t100\t.\tG\tA\t30\tPASS\t.\tGT\t0/1",
            "chr1\t200\t.\tC\tT\t30\tPASS\t.\tGT\t1/1",
        ],
    )
    assert list(read_biallelic_snvs(vcf_path)) == [
        SnvSite("chr1", 100, "G", "A"),
        SnvSite("chr1", 200, "C", "T"),
    ]
