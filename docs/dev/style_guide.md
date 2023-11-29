RNAvigate style guide
---------------------

RNAvigate adheres to the Google Python Style Guide found
[here](https://google.github.io/styleguide/pyguide.html). What follows is a
largely abridged version with some of my own additions regarding this project.

In general, keep it as simple as possible, human readable, and consistent.

TOC:
- [Commit messages](#commit-messages)
- [Linting](#linting)
- [Type Annotated Code](#type-annotated-code)
- [Long lines](#long-lines)
  - [long imports](#long-imports)
  - [long function calls](#long-function-calls)
- [Naming](#naming)
- [Main](#main)

## Commit messages

Commit messages are generally limited to 50 characters, which is very short.
These short-hand codes help and should be used in all commits. These codes
should be in parentheses, and parentheses should not be used elsewhere.

(+) Means that a feature was added.
(-) Means a feature or name was removed.
(~) Means a behavior was altered.
(fix) Means a bug was fixed
(doc) Means only changes are to documentation

For large changes, a longer message is needed. I don't use a specific format
for this. It should be concise and clear, and describe how and why the changes
were made.

## Linting

Run pylint over your code using the .pylintrc in the top-level directory. Fix
these issues regularly.

## Type Annotated Code

Most of RNAvigate is not type-hinted, but they are fine to use. I think these
are most readable when they are one argument per line with the closing
parenthesis also on it's own line. similar to
[long function calls](#long-function-calls).

# Long lines

Maximum line length is 80 characters, except:
- long imports, to keep them explicit
- long URLs, these should get their own line

## Long imports

In general, only modules and submodules should be imported. This usually means
that imports will fit on one line. In cases where multiple submodules are
imported, use grouping to wrap the line. You'll see this applied often in
RNAvigate's init files, where many names are imported explicitly.

```python
from a_really_long_module_name_with_multiple_submodules import (
    submodule1,
    submodule2,
    submodule3,
    )
```

## Long function calls

Function calls occasionally wrap to multiple lines to avoid the 80 character
limit. I prefer this layout with a closing parenthesis indented on it's own
line. Function calls under 80 characters are fine on one line.

Do this:
```python
short_function(short_arg='')
a_long_function_name_with_a_few_short_arguments(
    x=1, y=2, message="hello there",
    )
a_long_function_name_with_long_arguments(
    renormalization_method="percentile",
    renormalization_bounds=(0.95, 0.99),
    )
```

Not this: (weird indentation)
```python
short_function(
    short_arg='')
a_long_function_name_with_a_few_short_arguments(x=1, y=2,
                                                message="hello there")
a_long_function_name_with_long_arguments(renormalization_method="percentile",
                                         renormalization_bounds=(0.95, 0.99))
```

## Naming

- module_name
- package_name
- ClassName
- method_name
- ExceptionName
- function_name
- GLOBAL_CONSTANT_NAME
- global_var_name
- instance_var_name
- function_parameter_name
- local_var_name
- query_proper_noun_for_thing
- send_acronym_via_https.

avoid abbreviations. Except for:
- nt for nucleotide
- pos for position (1-indexed inclusive)
- idx for index (0-indexed, exclusive right)
- i and j for interactions or base-pair indices or nested enumerations
- in simple list comprehensions
- e for exception identifiers (except Error as e: Raise e)
- f for open file handlers in with statements
- mathematical notation in small scopes (r for radius, etc.)
  - spell out greek names (delta_t, not dt or d_t, etc.)

Always use a .py filename extension. Never use dashes.

don't include the type in variable names (id_to_name, not id_to_name_dict)

## Main

For now, RNAvigate is explicitly a tool for Jupyter and scripting, not command
line use, no `if __name__ == '__main__':` statements are necessary.
