import argparse
import re

def is_comment_line(line):
    line = line.strip()
    if len(line) == 0:
        return False
    else: return line[0] == '#'


def test_is_comment_line():
    assert is_comment_line("# This is a comment.")
    assert is_comment_line("    # This is a comment.")
    assert not is_comment_line("This is not  comment.")
    assert not is_comment_line("        This is not  comment.")


def collect_next_lines_with_logic(lines, logic):
    if len(lines) == 0:
        return []
    if not logic(lines[0]):
        return []
    return [lines[0]] + collect_next_lines_with_logic(lines[1:], logic)


def collect_next_comment_lines(lines):
    return collect_next_lines_with_logic(lines, lambda line: is_comment_line(line))


def collect_next_code_lines(lines):
    return collect_next_lines_with_logic(lines, lambda line: not is_comment_line(line))


def linearize_comment_lines(lines):
    newlines = [re.sub(r"^\s*# ", "", line).rstrip() + " " for line in lines]
    return "# " + "".join(newlines).strip()


def join_comment_lines(lines):

    def join_comment_lines_helper(lines):
        if len(lines) == 0:
            return []
        if lines[0].strip() == '':
            return [''] + join_comment_lines_helper(lines[1:])
        if is_comment_line(lines[0]):
            collected = collect_next_comment_lines(lines)
            k = len(collected)
            s = linearize_comment_lines(collected)
            return [s] + join_comment_lines_helper(lines[k:])
        else:
            collected = collect_next_code_lines(lines)
            k = len(collected)
            return collected + join_comment_lines_helper(lines[k:])

    header = collect_next_comment_lines(lines)  # the list of comments composing the header of the file
    k = len(header)

    return header + join_comment_lines_helper(lines[k:])


def wrap_comment_line(line, width):
    parts = list(filter(lambda x: x != None, re.split(r"(\[.*?\]\(.*?\))|(\`.*?\`)|(\[.*?\]:.*)|(\s+)", line)))
    lines = []
    s = parts[0]
    for part in parts[1:]:
        if s != '' and s[-1] == ' ' and part == ' ':
            continue
        if len(part) == 1:
            s += part
            continue
        if len(s) + len(part) <= width:
            s += part
            continue
        lines.append(s.strip())
        s = '# ' + part.strip()
    if s != '# ':
        lines.append(s.strip())
    return lines


def test_wrap_comment_line():
    line = "# In this link [Database](https://reaktoro.org/cpp/classReaktoro_1_1Database.html) and also here [Element](https://reaktoro.org/cpp/classReaktoro_1_1Element.html), you `can enter code` using ticks."
    lines = wrap_comment_line(line, 60)
    assert lines == [
        "# In this link",
        "# [Database](https://reaktoro.org/cpp/classReaktoro_1_1Database.html)",
        "# and also here",
        "# [Element](https://reaktoro.org/cpp/classReaktoro_1_1Element.html),",
        "# you `can enter code` using ticks."
    ]


def wrap(inputfile, outputfile, width):
    lines = [line.rstrip('\n') for line in open(inputfile)]
    lines = join_comment_lines(lines)
    newlines = []
    for line in lines:
        if is_comment_line(line):
            newlines += wrap_comment_line(line, width)
        else:
            newlines.append(line)

    with open(outputfile, "w") as f:
        for line in newlines:
            f.write(line + "\n")

def main():

    parser = argparse.ArgumentParser(
        description='A comment line wrapper for python files.')

    parser.add_argument('input',
        help='the python file to have comment lines formatted')

    parser.add_argument('--output',
        help='the formatted output python file')

    parser.add_argument('--width',
        type=int,
        default=121,
        help='the maximum line width of comment lines')

    args = parser.parse_args()

    if args.output == None:
        args.output = args.input

    wrap(args.input, args.output, args.width)


if __name__ == "__main__":
    main()
