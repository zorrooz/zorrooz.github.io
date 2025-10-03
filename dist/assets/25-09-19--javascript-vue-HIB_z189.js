const r=`---\r
title: "Vue.js 前端开发入门"\r
date: "2025-09-19"\r
author: "zorrooz"\r
tags: ["JavaScript", "Vue.js", "前端开发", "组件化", "响应式"]\r
draft: false\r
description: "Vue.js 框架在科学数据可视化前端开发中的应用"\r
---\r
\r
# Vue.js 前端开发入门\r
\r
## 基本语法\r
\r
\`\`\`javascript\r
// Vue 组件示例\r
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
## 组件通信\r
\r
- Props：父组件向子组件传递数据\r
- Events：子组件向父组件发送消息\r
- Vuex：全局状态管理`;export{r as default};
