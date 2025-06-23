"""
Module for stuff that is not specifically part of the feature-based analysis library
but is useful for scripts.
"""
import typing
__all__ = [
    "read_fasta",
    "read_regions_csv",
    "iter_nested"
]
def read_fasta(path):
    """Read a fasta file into a list of (header, sequence) tuples"""
    import warnings
    fasta_list = []
    with open(path, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:]
                fasta_list.append([header, []])
            elif fasta_list:
                fasta_list[-1][1].append(line)
            elif line:
                warnings.warn("Ignoring non-empty line before first fasta header: {}".format(line))
        for index, (header, seqlines) in enumerate(fasta_list):
            fasta_list[index] = (header, "".join(seqlines))
        return fasta_list
    
def read_regions_csv(path):
    """
    Read a csv file of the following form:
    ```
    ProteinID,RegionID,Start,Stop
    ...
    ```
    Where `ProteinID` and `RegionID` form a unique primary key and (`Start`, `Stop`)
    are 0-indexed coordinates of the IDR in the whole protein sequence.

    Output
    ------
    A doubly-nested dict of the form
    {`<ProteinID>`: {`<RegionID>`: (`<Start>`, `<Stop>`)}}
    """
    import csv

    return_value = {}
    with open(path, "r") as file:
        reader = csv.DictReader(file)
        fieldnames = reader.fieldnames
        if fieldnames is None:
            return return_value
        if len(fieldnames) != 4:
            raise ValueError("Expected csv to be of the form `ProteinID,RegionID,Start,Stop`, got these columns: {}".format(fieldnames))
        protid_col, regionid_col, start_col, stop_col = fieldnames
        for row in reader:
            protid = row[protid_col]
            regionid = row[regionid_col]
            start = row[start_col]
            try:
                start = int(start)
            except ValueError:
                raise ValueError("Expected third column of csv to be a numeric bound, got: {}".format(start))
            stop = row[stop_col]
            try:
                stop = int(stop)
            except ValueError:
                raise ValueError("Expected third column of csv to be a numeric bound, got: {}".format(start))
            if (entry := return_value.get(protid)) is None:
                entry = return_value[protid] = {}
            entry[regionid] = (start, stop)
    return return_value

def read_nested_csv(path, pkey_depth, *, group_multiple=False):
    """
    Read a csv file of the following form:
    ```
    Pkey1,Pkey2,...,PkeyN,Col1,Col2,...
    ...
    ```
    Where all `Pkey`s form a unique primary key.

    Arguments
    ---------

    `path` : path to csv file

    `pkey_depth` : number of primary keys

    `group_multiple` : if `True`, will use lists to store rows,
                        and multiple rows with the same pkeys will
                        be in the same list
    Output
    ------
    If `group_multiple=False`, a n-depth nested dict of the form
    {`<Pkey1>`: {`<Pkey2>`: ... {"Col1": `<Col1>`, "Col2": `<Col2>`, ...} ...}}

    If `group_multiple=False`, a n-depth nested dict of the form
    {`<Pkey1>`: {`<Pkey2>`: ... [{"Col1": `<Col1>`, "Col2": `<Col2>`, ...}, ...] ...}}
    """
    import csv

    return_value = {}
    with open(path, "r") as file:
        reader = csv.DictReader(file)
        fieldnames = reader.fieldnames
        if fieldnames is None:
            return return_value
        if len(fieldnames) < pkey_depth:
            raise ValueError("Expected at least {} columns, got these columns: {}".format(pkey_depth, fieldnames))
        *pkey_columns, last_pkey_column = fieldnames[:pkey_depth]
        for row in reader:
            entry = return_value
            for pkey_col in pkey_columns:
                pkey = row.pop(pkey_col)
                if (entry_next := entry.get(pkey)) is None:
                    entry_next = entry[pkey] = {}
                entry = entry_next
            pkey = row.pop(last_pkey_column)
            if group_multiple:
                if (entry_last := entry.get(pkey)) is None:
                    entry_last = entry[pkey] = []
                entry_last.append(row)
            else:
                entry[pkey] = row
    return return_value

def iter_nested(nested_dict, depth) -> typing.Generator[typing.Tuple[typing.Any, ...]]:
    """
    Iterate over a nested dictionary.

    `nested_dict` : A nested dictionary

    `depth` : A number at least one which represents how many
                levels of keys there are in the dict.

    Examples
    --------
    Given the value:
    ```
    value = {
        "A1": {"A2": 30, "B2": 40},
        "B1": {"A2": 14, "B2": 90},
        "C1": {"CORMORANT": 1023, "GOOSE": {"strange thing": "spooky"}}
    }
    ```
    the expression:
    ```
    for k1, k2, v in iter_nested(value, 2):
        print(k1, k2, v)
    ```
    yields:
    ```
    A1 A2 30
    A1 B2 40
    B1 A2 14
    B1 B2 90
    C1 CORMORANT 1023
    C1 GOOSE {"strange thing": "spooky"}
    ```
    And the expression:
    ```
    for key, v in iter_nested(value, 1):
        print(key, v)
    ```
    is equivalent to:
    ```
    for key, v in value.items():
        print(key, v)
    ```
    And hopefully this is clear how this extends to triply or quadruply-nested dictionaries.
    """
    if depth == 1:
        for key, value in nested_dict.items():
            yield key, value
    stack = []
    top_iter = iter(nested_dict.items())
    while True:
        try:
            nxt = next(top_iter)
        except StopIteration:
            if not stack:
                return
            top_iter, _ = stack.pop()
            depth += 1
            continue
        key, value = nxt
        if depth > 0:
            stack.append((top_iter, key))
            top_iter = iter(value.items())
            depth -= 1
            continue
        yield tuple([k for _, k in stack] + [key, value])
