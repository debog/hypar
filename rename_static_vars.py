#!/usr/bin/env python3
import os
import re

def rename_static_variables(filepath):
    """Rename static variables in a file to add s_ prefix."""

    try:
        with open(filepath, 'r') as f:
            content = f.read()
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return False

    original_content = content
    lines = content.split('\n')

    # Find static variable declarations (not functions)
    # Pattern: static type varname = value; or static type varname;
    # Exclude function declarations (those with parentheses after the name)
    static_var_pattern = r'(\s*static\s+(?:const\s+)?(?:unsigned\s+)?(?:int|double|float|char|long|short|size_t)\s+)([a-zA-Z_][a-zA-Z0-9_]*)(\s*(?:=|;))'

    # Keep track of renamed variables
    renamed_vars = {}

    # First pass: find all static variables and determine new names
    for i, line in enumerate(lines):
        # Skip if it's a function declaration (contains '(' after variable name)
        if '(' in line and ')' in line:
            continue

        match = re.search(static_var_pattern, line)
        if match:
            var_name = match.group(2)
            # Skip if already has s_ prefix
            if not var_name.startswith('s_'):
                new_name = 's_' + var_name
                renamed_vars[var_name] = new_name
                # Replace in the declaration line
                lines[i] = lines[i].replace(match.group(2), new_name, 1)

    if not renamed_vars:
        return False

    # Second pass: replace all usages of renamed variables
    # Be careful to only replace whole words
    for old_name, new_name in renamed_vars.items():
        # Use word boundaries to avoid partial matches
        pattern = r'\b' + re.escape(old_name) + r'\b'
        content = '\n'.join(lines)
        content = re.sub(pattern, new_name, content)
        lines = content.split('\n')

    # Write back
    new_content = '\n'.join(lines)
    if new_content != original_content:
        try:
            with open(filepath, 'w') as f:
                f.write(new_content)
            print(f"Renamed {len(renamed_vars)} variable(s) in: {filepath}")
            for old, new in renamed_vars.items():
                print(f"  {old} -> {new}")
            return True
        except Exception as e:
            print(f"Error writing {filepath}: {e}")
            return False

    return False

def main():
    src_dir = '/home/ghosh/Codes/hypar/src'
    include_dir = '/home/ghosh/Codes/hypar/include'

    modified_count = 0

    # Process source files
    for root, dirs, files in os.walk(src_dir):
        for filename in files:
            if filename.endswith('.c') or filename.endswith('.cpp') or filename.endswith('.cu'):
                filepath = os.path.join(root, filename)
                if rename_static_variables(filepath):
                    modified_count += 1

    # Process header files
    for root, dirs, files in os.walk(include_dir):
        for filename in files:
            if filename.endswith('.h') or filename.endswith('.hpp'):
                filepath = os.path.join(root, filename)
                if rename_static_variables(filepath):
                    modified_count += 1

    print(f"\n{'='*60}")
    print(f"Total files modified: {modified_count}")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()
