# Usage
Demo with `python plot.py -f -i example -o example`.

Call with
`python plot.py -f -i filename -o out -s 1`,
where `filename` (optional) is the prefix for the input files described below.
If `filename` is not provided, the prefix will default to `input`.
See the `example.<xy,nodelist,edgelist>` files for how to format the input files.

### Options
- `-f --full`:
    - boolean flag
    - when set, makes a full, stand-alone document
    - otherwise, makes only the content between `\begin{tikzpicture}` `\end{tikzpicture}` tags
- `-n --nodecolor`:
    - boolean flag
    - when set, colors edges according to labels in the last column of the `filename.nodelist` input
- `-e --edgecolor`:
    - boolean flag
    - when set, colors edges according to labels in the last column of the `filename.edgelist` input
- `-i --input`:
    - specifies the prefix for the input filename
    - defaults to `input`
- `-o --output`:
    - specifies the prefix for the output filename
    - defaults to `output`
- `-s --scale`:
    - scaling factor on the inter-node separation of the resultant tikz picture
    - defaults to `1`

## Input file formats
- `filename.xy`:
    - space-delimited file where line `i` contains two real numbers `xᵢ yᵢ`
    - `xᵢ yᵢ` are the Cartesian coordinates of the node on line `i` in `filename.nodelist`
- `filename.nodelist`:
    - file where line `i` contains a string `vᵢ` and (optionally) a label `l ∈ {-1, 0, 1, 2, 3, 4, 5, 6}`
    - the string `vᵢ` is the name of the node as used in the `filename.edgelist` file
    - the label `l` is used for node coloring, if enabled
- `filename.edgelist`:
    - space-delimited file where line `i` two node names `vᵢ vⱼ` and (optionally) a label `l ∈ {-1, 0, 1}`
    - the node names `vᵢ` and `vⱼ` should correspond to the names given in `filename.nodelist`
    - the label `l` is used for edge coloring, if enabled

## Output
A LaTeX file named `output.tex`, where `output` is the value set by the `-o` flag, containing a tikz picture of the graph defined by the input files.
