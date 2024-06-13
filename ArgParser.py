import argparse



def getArgParser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-f",
        "--input-file",
        type=str,
        help="Input file [.fna]",
        required=True
    )
    parser.add_argument(
        "-n",
        "--nucleotide-quality",
        type=float,
        help="Average nucleotide quality",
        required=True
    )
    parser.add_argument(
        "-s",
        "--nucleotide-stddev",
        type=float,
        help="Average nucleotide quality",
        required=True
    )
    parser.add_argument(
        "-c",
        "--coverage",
        type=float,
        help="nucleotide coverage",
        required=True,
    )
    parser.add_argument(
        "-m",
        "--snv-error-rate",
        type=float,
        help="Single Nucleotide Variant mutation error rate",
        default=0.0
    )
    parser.add_argument(
        "-i",
        "--insertion-error-rate",
        type=float,
        help="Insertion error rate",
        default=0.0
    )
    parser.add_argument(
        "-d",
        "--deletion-error-rate",
        type=float,
        help="Deletion error rate",
        default=0.0
    )
    parser.add_argument(
        "-o",
        "--output-sam-file",
        type=str,
        help="Output sam file", #path (output .sam file)
        required=True
    )
    parser.add_argument(
        "--output-fq-file1",
        type=str,
        help="Output first fq file path",  # path (output .sam file)
        required=True
    )
    parser.add_argument(
        "--output-fq-file2",
        type=str,
        help="Output second fq file path",  # path (output .sam file)
        required=True
    )
    return parser

if __name__ == '__main__':
    argParser = getArgParser()
    args = argParser.parse_args()
    print(args)