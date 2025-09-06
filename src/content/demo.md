# Markdown内容demo

这是一段普通文本，**加粗文本**，*斜体文本*，以及[链接](https://vuejs.org)。

## 代码高亮

行内代码是：`console.log('Hello, Vue Markdown!');`。

```javascrip
// JavaScript代码示例
function greet(name) {
  return 'Hello, ' + name + '!';
}

console.log(greet('Vue Markdown'));
```

## 数学公式

行内公式：$E = mc^2$

块级公式：
$$E=mc^2$$

$$f'(x) = \lim_{h \to 0} \frac{f(x+h) - f(x)}{h}$$

$$
\int_a^b f(x) \, \mathrm{d}x = F(b) - F(a)
$$

## 表格示例

| 功能        | 支持情况   | 备注          |
| ----------- | ---------- | ------------- |
| 标题        | 支持       | 最多支持六级标题 |
| 代码块      | 支持       | 支持语法高亮    |
| 数学公式    | 支持       | 使用KaTeX渲染  |
| 表格        | 支持       | 标准表格语法    |
| 图片        | 支持       | 本地和网络图片  |

## 图片示例

![Vue Logo](https://vuejs.org/images/logo.png)

## 引用块
>
> 这是一个引用块示例
> 它可以包含多行内容

## 列表

### 有序列表

1. 第一项
2. 第二项
3. 第三项

### 无序列表

- Vue
- Markdown
- KaTeX
- Highlight.js

```python
def greet(name):
    return 'Hello, ' + name + '!'

print(greet('Vue Markdown'))
1 + 1
```

---
以上是Markdown渲染的基本功能演示。
