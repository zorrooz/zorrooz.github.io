const n=`# Markdown 语法示例

以下是对常见 Markdown 语法功能的演示。

## 1. 基础文本样式

### 三级标题

#### 四级标题

##### 五级标题

###### 六级标题

普通段落文本示例。

1. 粗体文本
    - **星号形式**
    - **下划线形式**
2. 斜体文本
    - *星号形式*
    - *下划线形式*

## 2. 链接与图片

这是一个行内链接示例：[Vue 官网](https://vuejs.org)

以下是嵌入图片的示例：  
![Vue Logo](https://vuejs.org/images/logo.png)

## 3. 块引用

> 这是一个块引用。  
> 块引用常用于强调某些内容或展示他人的陈述。

## 4. 表格与行内代码

| 样式类型          | 示例 |
| ----------------- | ---- |
| 行内代码          | \`console.log("Hello World")\` |
| 行内数学公式      | $f'(x) = \\lim_{h \\to 0} \\frac{f(x+h) - f(x)}{h}$ |

## 5. 代码块

缩进式代码块：

    console.log("Hello_World");

围栏式代码块：

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

## 6. 数学公式

块级数学公式示例：

$$
\\int_a^b f(x) \\, \\mathrm{d}x = F(b) - F(a)
$$

---
以上是常用的 Markdown 语法功能示例。
`;export{n as default};
