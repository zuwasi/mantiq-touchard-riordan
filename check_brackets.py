import sys
path = r"C:\Amp_demos\Eli-Bango-JCT\mantiq-project\src\SetPartitionQSym.wl"
text = open(path, encoding="utf-8").read()
pairs = {"[": "]", "(": ")", "{": "}"}
closes = {"]": "[", ")": "(", "}": "{"}
stack = []
line = 1
i = 0
n = len(text)
while i < n:
    ch = text[i]
    if ch == "\n":
        line += 1
        i += 1
        continue
    if text.startswith("(*", i):
        depth = 1
        i += 2
        while i < n and depth:
            if text.startswith("(*", i):
                depth += 1; i += 2
            elif text.startswith("*)", i):
                depth -= 1; i += 2
            else:
                if text[i] == "\n": line += 1
                i += 1
        continue
    if ch == '"':
        i += 1
        while i < n and text[i] != '"':
            if text[i] == "\n": line += 1
            i += 1
        i += 1
        continue
    if ch in pairs:
        stack.append((pairs[ch], line))
    elif ch in closes:
        if not stack:
            print("extra close at line", line, ch); sys.exit()
        top, ln = stack.pop()
        if top != ch:
            print("mismatch at line", line, "expected", top, "got", ch, "opened at", ln); sys.exit()
    i += 1
print("unclosed:" if stack else "balanced", stack[-12:] if stack else "")
