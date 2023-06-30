import sys

if __name__ == "__main__":

    n = int(sys.argv[1])

    f = open("clique_%d.txt" % n, "w")
    for i in range(0, n):
        for j in range(i + 1, n):
            f.write("%d %d" % (i, j))
            if (i < (n - 1) or j < (n - 1)):
                f.write("\n")
    f.close()
