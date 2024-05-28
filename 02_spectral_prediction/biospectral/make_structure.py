import sys

from biospectral.generate_structures import main

if __name__ == "__main__":
    try:
        main()

    except KeyboardInterrupt:
        print("Interrupted")
        sys.exit(0)
