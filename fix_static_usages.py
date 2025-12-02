#!/usr/bin/env python3
import os
import re

def fix_static_variable_usages(filepath):
    """Fix usages of static variables in files where they're declared."""

    try:
        with open(filepath, 'r') as f:
            content = f.read()
    except:
        return False

    original_content = content

    # Find static variables declared with s_ prefix in this file
    # Pattern: static ... s_varname = ...
    static_decl_pattern = r'static\s+(?:const\s+)?(?:\w+\s+)+s_(\w+)\s*='

    declared_vars = set()
    for match in re.finditer(static_decl_pattern, content):
        var_name = match.group(1)
        declared_vars.add(var_name)

    if not declared_vars:
        return False

    # Now replace usages of these variables (without s_ prefix)
    for var_name in declared_vars:
        # Match the variable without s_ prefix, using word boundaries
        # But don't match if it's already prefixed with s_
        old_pattern = r'(?<!s_)\b' + re.escape(var_name) + r'\b'
        new_name = 's_' + var_name
        content = re.sub(old_pattern, new_name, content)

    if content != original_content:
        try:
            with open(filepath, 'w') as f:
                f.write(content)
            print(f"Fixed usages in: {filepath}")
            print(f"  Variables: {', '.join(sorted(declared_vars))}")
            return True
        except:
            return False

    return False

def main():
    src_dir = '/home/ghosh/Codes/hypar/src'
    include_dir = '/home/ghosh/Codes/hypar/include'

    modified_count = 0

    # Process source files
    for root, dirs, files in os.walk(src_dir):
        for filename in files:
            if filename.endswith(('.c', '.cpp', '.cu')):
                filepath = os.path.join(root, filename)
                if fix_static_variable_usages(filepath):
                    modified_count += 1

    # Process header files
    for root, dirs, files in os.walk(include_dir):
        for filename in files:
            if filename.endswith(('.h', '.hpp')):
                filepath = os.path.join(root, filename)
                if fix_static_variable_usages(filepath):
                    modified_count += 1

    print(f"\n{'='*60}")
    print(f"Total files modified: {modified_count}")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()
