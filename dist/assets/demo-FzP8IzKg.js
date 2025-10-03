const r=`# Markdown 语法示例\r
\r
以下是对常见 Markdown 语法功能的演示。\r
\r
## 1. 基础文本样式\r
\r
### 三级标题\r
\r
#### 四级标题\r
\r
##### 五级标题\r
\r
###### 六级标题\r
\r
普通段落文本示例。\r
\r
1. 粗体文本\r
    - **星号形式**\r
    - **下划线形式**\r
2. 斜体文本\r
    - *星号形式*\r
    - *下划线形式*\r
\r
## 2. 链接与图片\r
\r
这是一个行内链接示例：[Vue 官网](https://vuejs.org)\r
\r
以下是嵌入图片的示例：  \r
![Vue Logo](https://vuejs.org/images/logo.png)\r
\r
## 3. 块引用\r
\r
> 这是一个块引用。  \r
> 块引用常用于强调某些内容或展示他人的陈述。\r
\r
## 4. 表格与行内代码\r
\r
| 样式类型          | 示例 |\r
| ----------------- | ---- |\r
| 行内代码          | \`console.log("Hello World")\` |\r
| 行内数学公式      | $f'(x) = \\lim_{h \\to 0} \\frac{f(x+h) - f(x)}{h}$ |\r
\r
## 5. 代码块\r
\r
缩进式代码块：\r
\r
    console.log("Hello_World");\r
\r
围栏式代码块：\r
\r
\`\`\`html\r
<style>\r
.markdown-body img {\r
  max-width: 100%;\r
  height: auto;\r
  display: block;\r
  margin: 1em auto;\r
}\r
</style>\r
\`\`\`\r
\r
## 6. 数学公式\r
\r
块级数学公式示例：\r
\r
$$\r
\\int_a^b f(x) \\, \\mathrm{d}x = F(b) - F(a)\r
$$\r
\r
---\r
以上是常用的 Markdown 语法功能示例。\r
`;export{r as default};
