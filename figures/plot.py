import argparse

def plot(full: bool, in_filename: str, out_filename: str, scale: float = 1):
    with open(f'{in_filename}.xy', 'r') as cfile, \
            open(f'{in_filename}.nodelist', 'r') as nfile, \
            open(f'{in_filename}.edgelist', 'r') as efile:
        in_xy = [tuple(map(float, line.strip().split(' '))) for line in cfile]
        in_nodes = [line.strip() for line in nfile]
        in_edges = [tuple(line.strip().split(' ')) for line in efile]

    indent = '    ' if full else ''

    predoc = (
        '\\documentclass{standalone}\n' +
        '\\include{tikz.tex}\n' +
        '\\begin{document}\n'
    ) if full else ''
    pre = indent + '\\begin{tikzpicture}' + f'[scale={scale}]\n'
    nodes = ''.join(f'''{indent}    \\node [node] at ({x}, {y}) ({v}) {{}};\n'''
                    for (x, y), v in zip(in_xy, in_nodes))
    edges = ''.join(f'''{indent}    \\draw [edge] ({u}) to ({v});\n'''
                    for u, v in in_edges)
    post = indent + '\\end{tikzpicture}'
    postdoc = '\n\\end{document}' if full else ''

    with open(f'{out_filename}.tex', 'w') as outfile:
        outfile.write(predoc + pre + nodes + edges + post + postdoc)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f',
                      '--full',
                      action='store_true',
                      help = 'Whether to make a full, stand-alone document')
    parser.add_argument('-i',
                      '--input',
                      type=str,
                      default='input',
                      help = 'Input file prefix')
    parser.add_argument('-o',
                      '--output',
                      type=str,
                      default='output',
                      help = 'Output file prefix')
    parser.add_argument('-s',
                      '--scale',
                      type=float,
                      default=1,
                      help = 'Scale for the tikz figure')
    args = parser.parse_args()

    plot(args.full, args.input, args.output, scale=args.scale)
