"""
An interface for including any python function with as a function in the [`Featurizer`].

Whereas a native or builtin feature can be specified in json via something like:
```
{
    "my-count-feature": {
        "compute": "count",
        "pattern": "[ILVM]"
    },
    "my-lr-feature": {
        "compute": "log_ratio",
        "numerator": "N",
        "denominator": "M"
    }
}

You can specify a custom feature by setting `compute=custom`. The required information in
the rest of the dict is given in the docs for [`compile_custom_featurizer`].
```
"""
def compile_custom_featurizer(features_dict):
    """
    A utility function for specifying custom features in other python scripts.

    Schema
    ------
    `features_dict` should be a dict of feature names to json objects with three
    fields:
    1. `libpath`, the path where the script is
    2. `funcname`, the name of the function in the globals of the script
    3. `kwargs`, optional keyword arguments to curry into the function

    Example
    -------
    Here's an example feature dict:
    ```
    {
        "my-script-feature": {
            "libpath": "my-script-path.py",
            "funcname": "compute",
            "compute": "custom",
            "kwargs": {
                "kw_a": 5
            }
        }
    }
    ```
    with the corresponding python file at `my-script-path.py`
    ```
    def compute(sequence: str, kw_a: int):
        return len(sequence) + kw_a
    ```
    """
    import os
    import runpy
    import functools

    return_value = {}
    runpy_results = {}
    errors = {}
    for featname, feature_params in features_dict.items():
        if (libpath := feature_params.get("libpath")) is None:
            errors[featname] = ValueError(
                "feature needs a script path as a `libpath` field"
            )
            continue
        if not isinstance(libpath, str):
            errors[featname] = TypeError("expected a script path string at `libpath`")
            continue
        if (funcname := feature_params.get("funcname")) is None:
            errors[featname] = ValueError(
                "feature needs a function name as a `funcname` field"
            )
            continue
        if not isinstance(funcname, str):
            errors[featname] = ValueError(
                "expected a function name string as a `funcname` field"
            )
            continue
        keypath = os.path.realpath(libpath)
        if (module_or_error := runpy_results.get(keypath)) is None:
            try:
                module_or_error = runpy_results[keypath] = runpy.run_path(libpath)
            except (ImportError, OSError) as e:
                errors[featname] = runpy_results[keypath] = e
                continue
        if isinstance(module_or_error, Exception):
            errors[featname] = module_or_error
            continue
        module = module_or_error
        if (func := module.get(funcname)) is None:
            errors[featname] = RuntimeError(
                "No function `%s` in module at `%s`" % (funcname, libpath)
            )
            continue
        if (kwargs := feature_params.get("kwargs")) is None:
            kwargs = {}
        return_value[featname] = functools.partial(func, **kwargs)
    return return_value, errors
