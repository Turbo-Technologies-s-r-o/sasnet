"""Command line: python -m sasnet network.json [report.md]"""
import sys

from .network import Network
from .report import markdown


def main(argv=None):
    argv = argv or sys.argv[1:]
    if not argv:
        print(__doc__)
        return 2
    net = Network.from_json(argv[0])
    r = net.solve()
    text = markdown(net, r, title=argv[0])
    if len(argv) > 1:
        with open(argv[1], "w", encoding="utf-8") as fh:
            fh.write(text)
    print(text)
    return 0 if r.converged else 1


if __name__ == "__main__":
    sys.exit(main())
