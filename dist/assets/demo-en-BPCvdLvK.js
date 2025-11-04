const n=`# Markdown Syntax Examples

Below is a demonstration of common Markdown syntax features.

## 1. Basic Text Styles

### Level 3 Heading

#### Level 4 Heading

##### Level 5 Heading

###### Level 6 Heading

Regular paragraph text example.

1. Bold Text
    - **Asterisk format**
    - **Underscore format**
2. Italic Text
    - *Asterisk format*
    - *Underscore format*

## 2. Links and Images

This is an inline link example: [Vue Official Website](https://vuejs.org)

Below are embedded image examples:  
![Vue Logo](https://vuejs.org/images/logo.png)

## 3. Block Quotes

> This is a block quote.  
> Block quotes are commonly used to emphasize certain content or display others' statements.

## 4. Tables and Inline Code

| Style Type          | Example |
| ------------------- | ------- |
| Inline Code         | \`console.log("Hello World")\` |
| Inline Math Formula | $f'(x) = \\lim_{h \\to 0} \\frac{f(x+h) - f(x)}{h}$ |

## 5. Code Blocks

Indented code block:

    console.log("Hello_World");

Fenced code block:

\`\`\`html
<style>
.markdown-body img {
  max-width: 100%;
  height: auto;
  display: block;
  margin: 1em auto;
}
</style>
\`\`\`

## 6. Mathematical Formulas

Block-level mathematical formula example:

$$
\\int_a^b f(x) \\, \\mathrm{d}x = F(b) - F(a)
$$

---
Above are examples of commonly used Markdown syntax features.`;export{n as default};
