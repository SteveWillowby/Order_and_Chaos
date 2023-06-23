import sys

if __name__ == "__main__":

    if (len(sys.argv) != 3):
        print("Error! Takes exactly 2 arguments.")
        exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    f1 = open(input_file, "r")
    edge_list = f1.read()
    f1.close()

    f2 = open(output_file, "w")
    f2.write(edge_list.replace(" ", ","))
    f2.close()
