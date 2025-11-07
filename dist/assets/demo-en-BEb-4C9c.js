const n=`# Markdown Syntax Examples

The following demonstrates common Markdown syntax features.

## 1. Basic Text Styles

### Level 3 Heading

#### Level 4 Heading

##### Level 5 Heading

###### Level 6 Heading

Example of regular paragraph text.

1. Bold Text
    - **Asterisk form**
    - **Underscore form**
2. Italic Text
    - *Asterisk form*
    - *Underscore form*

## 2. Links and Images

Here is an inline link example: [Vue Official Website](https://vuejs.org)

Below is an example of embedding an image:  
![Vue Logo](https://vuejs.org/images/logo.png)

## 3. Blockquotes

> This is a blockquote.  
> Blockquotes are often used to emphasize certain content or display someone else's statement.

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

Block-level math formula example:

$$
\\int_a^b f(x) \\, \\mathrm{d}x = F(b) - F(a)
$$

---
The above are examples of commonly used Markdown syntax features.`;export{n as default};
