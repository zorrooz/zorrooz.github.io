const n=`---\r
title: "Vue.js Frontend Development Introduction"\r
date: "2025-09-19"\r
author: "zorrooz"\r
tags: ["JavaScript", "Vue.js", "Frontend Development", "Componentization", "Reactive"]\r
draft: false\r
description: "Application of Vue.js framework in scientific data visualization frontend development"\r
---\r
\r
# Vue.js Frontend Development Introduction\r
\r
## Basic Syntax\r
\r
\`\`\`javascript\r
// Vue component example\r
export default {\r
  data() {\r
    return {\r
      message: 'Hello Vue!'\r
    }\r
  },\r
  methods: {\r
    updateMessage() {\r
      this.message = 'Updated!'\r
    }\r
  }\r
}\r
\`\`\`\r
\r
## Component Communication\r
\r
- Props: Parent component passes data to child component\r
- Events: Child component sends messages to parent component\r
- Vuex: Global state management\r
\r
Below is a demonstration of common Markdown syntax features.\r
\r
## 1. Basic Text Styles\r
\r
### Level 3 Heading\r
\r
#### Level 4 Heading\r
\r
##### Level 5 Heading\r
\r
###### Level 6 Heading\r
\r
Regular paragraph text example.\r
\r
1. Bold Text\r
    - **Asterisk format**\r
    - **Underscore format**\r
2. Italic Text\r
    - *Asterisk format*\r
    - *Underscore format*\r
\r
## 2. Links and Images\r
\r
This is an inline link example: [Vue Official Website](https://vuejs.org)\r
\r
Below are embedded image examples:\r
\r
 ![city](./city.png)\r
![Vue Logo](https://vuejs.org/images/logo.png)\r
\r
## 3. Block Quotes\r
\r
> This is a block quote.  \r
> Block quotes are commonly used to emphasize certain content or display others' statements.\r
\r
## 4. Tables and Inline Code\r
\r
| Style Type          | Example |\r
| ------------------- | ------- |\r
| Inline Code         | \`console.log("Hello World")\` |\r
| Inline Math Formula | $f'(x) = \\lim_{h \\to 0} \\frac{f(x+h) - f(x)}{h}$ |\r
\r
## 5. Code Blocks\r
\r
Indented code block:\r
\r
    console.log("Hello_World");\r
\r
Fenced code block:\r
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
## 6. Mathematical Formulas\r
\r
Block-level mathematical formula example:\r
\r
$$\r
\\int_a^b f(x) \\, \\mathrm{d}x = F(b) - F(a)\r
$$\r
\r
---\r
Above are examples of commonly used Markdown syntax features.`;export{n as default};
