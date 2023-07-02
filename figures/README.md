# usage
Demo with `python plot.py -f -i example -o example`.

Call with
`python plot.py -f -i in -o out -s 1`,
where `in` (optional) is the prefix for the input files described below.
If `fn` is not provided, the prefix will default to `input`.
See the `example.<xy,nodelist,edgelist>` files for how to format the input files.

### options
- `-f --full`:
    - boolean flag
    - when set, makes a full, stand-alone document
    - otherwise, makes only the content between `\begin{tikzpicture}` `\end{tikzpicture}` tags
- `-i --input`:
    - specifies the prefix for the input filename
    - defaults to `input`
- `-o --output`:
    - specifies the prefix for the output filename
    - defaults to `output`
- `-s --scale`:
    - scaling factor on the inter-node separation of the resultant tikz picture
    - defaults to `1`

## input file format
- `fn.xy`:
    - space-delimited file where line `i` contains two real numbers `xᵢ yᵢ`
    - `xᵢ yᵢ` are the Cartesian coordinates of the node on line `i` in `fn.nodelist`
- `fn.nodelist`:
    - file where line `i` contains the name of node `vᵢ` as a contiguous string
- `fn.edgelist`:
    - space-delimited file where line `i` contains the name of node `vᵢ`

## output
A latex file named `output.tex`, where `output` is the value set by the `-o` flag, containing a .
