""" Command line code """
import logging
import argparse
from pyhlamsa import HLAmsaEX, KIRmsa, CYPmsa, Genemsa, Familymsa

logger = logging.getLogger()
logging.basicConfig(level=logging.INFO)


def add_parser() -> argparse.ArgumentParser:
    """Add two command: download and view"""
    parser = argparse.ArgumentParser(
        prog="pyHLAMSA",
        description="Read hla/cyp/kir, do somthing and save to fasta/bam/vcf",
    )
    parser.add_argument("--debug", action="store_true", help="Debug Output")

    subparser = parser.add_subparsers(title="subcommand", dest="subcommand")
    # download command
    parser_download = subparser.add_parser(
        "download",
        help="Download hla/cyp/kir and save to msa",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_download.add_argument(
        "--family",
        required=True,
        choices=["hla", "kir", "cyp"],
        help="Type of gene family to download",
    )
    parser_download.add_argument(
        "--seq-type",
        choices=["gen", "nuc", "merged"],
        default="gen",
        help="The sequence type for parsing",
    )
    parser_download.add_argument(
        "--version", default="Latest", help="Version of database"
    )
    parser_download.add_argument(
        "--db-folder",
        default="",
        help="Folder of database "
        "If folder not exist, it will download "
        "the specific version databse to this folder",
    )
    parser_download.add_argument("name", help="Prefix of output name")
    parser_download.add_argument(
        "--include-genes",
        nargs="+",
        default=None,
        help="Specific the genes you want to parse",
    )
    parser_download.add_argument(
        "--consensus-name",
        default="*consensus",
        help="The name of the new consensus sequences "
        "automatically generated "
        "(format: {gene}{consensus-name})",
    )
    # view command
    parser_view = subparser.add_parser(
        "view",
        help="View the msa and save to file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_view.add_argument("input_name", help="Prefix name of .json and .fa")
    parser_view.add_argument(
        "--name", default=None, required=False, help="Prefix name of output format"
    )
    parser_view.add_argument(
        "--no-show", action="store_true", help="Don't show the segment of MSA to stdout"
    )
    parser_view.add_argument(
        "--show-diff",
        action="store_true",
        help="Only show the base has variantion in MSA",
    )
    parser_view.add_argument("--save", action="store_true", help="Save to .json and fa")
    parser_view.add_argument(
        "--region",
        nargs="+",
        help="Specific the region (split by space)" "Example: intron1 exon2",
    )
    parser_view.add_argument(
        "--position", type=str, help="Specific the position. Example: 100-400"
    )
    parser_view.add_argument(
        "--include-alleles", nargs="+", default=None, help="Show these alleles only"
    )
    parser_view.add_argument("--exclude-alleles", nargs="+", help="Remove some alleles")
    parser_view.add_argument(
        "--fasta-msa",
        action="store_true",
        help="Save to fasta (All sequences are equal length)",
    )
    parser_view.add_argument(
        "--fasta-gapless", action="store_true", help="Save to fasta"
    )
    parser_view.add_argument("--bam", action="store_true", help="Save to bamfile")
    parser_view.add_argument("--gff", action="store_true", help="Save to gff")
    parser_view.add_argument("--vcf", action="store_true", help="Save to vcf")
    return parser


def download_command(args: argparse.Namespace):
    """ " Download database and save to msa"""
    if args.seq_type == "merged":
        filetype = ["nuc", "gen"]
    else:
        filetype = [args.seq_type]

    if args.family.lower() == "hla":
        family: Familymsa = HLAmsaEX(
            genes=args.include_genes,
            version=args.version,
            filetype=filetype,
            imgt_folder=args.db_folder,
        )
    elif args.family.lower() == "kir":
        family = KIRmsa(
            genes=args.include_genes,
            version=args.version,
            filetype=filetype,
            ipd_folder=args.db_folder,
        )
    elif args.family == "cyp":
        family = CYPmsa(
            genes=args.include_genes,
            version=args.version,
            filetype=filetype,
            pharmvar_folder=args.db_folder,
        )
    for gene_name, msa in family.items():
        msa = msa.shrink()
        msa = msa.reset_index()
        msa.append(f"{gene_name}{args.consensus_name}", msa.get_consensus())
        msa.set_reference(f"{gene_name}{args.consensus_name}")
        msa.save_msa(f"{args.name}.{gene_name}.fa", f"{args.name}.{gene_name}.json")
        logger.info(f"Save to {args.name}.{gene_name}.*")


def extract_msa(args) -> Genemsa:
    """read msa, selection region/position/alleles"""
    msa = Genemsa.load_msa(f"{args.input_name}.fa", f"{args.input_name}.json")
    if args.include_alleles is not None:
        msa = msa.select_allele(args.include_alleles)
    if args.exclude_alleles:
        msa = msa.select_allele(
            list(set(msa.list_alleles()) - set(args.exclude_alleles))
        )
    if args.region:
        msa = msa.select_block(args.region)
    if args.position:
        # format: start-end
        # string -> int -> 1-base-position
        start, end = map(lambda i: i - 1, map(int, args.position.split("-")))
        msa = msa[start:end]
    return msa


def write_to_files(args, msa: Genemsa):
    """save msa to another format"""
    if not args.name:
        if (
            args.bam
            or args.gff
            or args.fasta_gapless
            or args.fasta_msa
            or args.vcf
            or args.save
        ):
            raise ValueError("No output name specific")
        return
    if args.input_name == args.name:
        raise ValueError(
            "output name should not same as input name," "it will overwrite it"
        )
    save_ref_seq = False
    if args.save:
        msa.save_msa(f"{args.name}.fa", f"{args.name}.json")
        logger.info(f"Save to {args.name}.fa {args.name}.json")
    msa = msa.shrink()
    if args.bam:
        save_ref_seq = True
        msa.to_bam(f"{args.name}.bam")
        logger.info(f"Save to {args.name}.bam")
    if args.gff:
        save_ref_seq = True
        msa.to_gff(f"{args.name}.gff")
        logger.info(f"Save to {args.name}.gff")
    if args.fasta_gapless:
        msa.to_fasta(f"{args.name}.nogap.fa", gap=False)
        logger.info(f"Save to {args.name}.nogap.fa")
    if args.fasta_msa:
        msa.to_fasta(f"{args.name}.msa.fa", gap=True)
        logger.info(f"Save to {args.name}.msa.fa")
    if args.vcf:
        save_ref_seq = True
        msa.to_vcf(f"{args.name}.vcf.gz")
        logger.info(f"Save to {args.name}.vcf.gz")
    if save_ref_seq:
        msa.to_fasta(f"{args.name}.ref.fa", gap=False, ref_only=True)
        logger.info(f"Save to {args.name}.ref.fa")


def run_command(args: argparse.Namespace):
    """Run command by args"""
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    # download databse
    if args.subcommand == "download":
        download_command(args)
    # read msa and transfer to fasta/bam/vcf/gff
    elif args.subcommand == "view":
        msa = extract_msa(args)
        if not args.no_show:
            if args.show_diff:
                msa.print_snv()
            else:
                msa.print_alignment_diff()
        write_to_files(args, msa)


def main():
    """Main function for command line"""
    parser = add_parser()
    args = parser.parse_args()
    logger.debug(args)
    if args.subcommand:
        run_command(args)
    else:
        print(parser.print_help())
