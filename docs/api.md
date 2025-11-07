````markdown
# API documentation

This page outlines how to generate API docs from the source using Sphinx autodoc.

1. Ensure your package is importable by Sphinx (the `conf.py` adds the repository root to sys.path).
2. Add autodoc directives to this file, for example:

```rst
.. automodule:: hypar
    :members:
    :undoc-members:
    :show-inheritance:
```

If you prefer Markdown style with MyST, you can embed an autodoc directive:

```markdown
```{automodule} hypar
:members:
:undoc-members:
:show-inheritance:
```
```

(When adding these directives in Markdown with MyST, use the fenced directive syntax as shown above.)
````