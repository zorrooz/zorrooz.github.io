<!-- Markdown.vue -->
<template>
  <div class="markdown-viewer">
    <!-- 目录 -->
    <div v-if="tocHtml" v-html="tocHtml" class="toc mb-4 p-3 bg-light rounded shadow-sm"></div>

    <!-- 内容 -->
    <div class="markdown-body" v-html="html"></div>
  </div>
</template>

<script>
import MarkdownIt from 'markdown-it';
import MarkdownItToc from 'markdown-it-toc-done-right';
import MarkdownItKatex from 'markdown-it-katex';
import hljs from 'highlight.js';

// 注意：CSS 在 main.js 或全局样式中引入（见下方说明）
// import 'katex/dist/katex.min.css'
// import 'highlight.js/styles/github.css'

export default {
  name: 'MarkdownViewer',
  props: {
    content: {
      type: String,
      required: true,
    },
  },
  data() {
    return {
      html: '',
      tocHtml: '',
    };
  },
  watch: {
    content: {
      immediate: true,
      handler() {
        this.render();
      },
    },
  },
  methods: {
    render() {
      // 初始化 markdown-it
      const md = new MarkdownIt({
        html: true,
        linkify: true,
        typographer: true,
        highlight: (str, lang) => {
          if (lang && hljs.getLanguage(lang)) {
            try {
              return (
                `<pre class="hljs"><code>` +
                hljs.highlight(str, { language: lang }).value +
                `</code></pre>`
              );
            } catch (error) {
              console.error('代码高亮错误:', error);
            }
          }
          // 默认回退：转义 HTML
          return (
            `<pre class="hljs"><code>` +
            md.utils.escapeHtml(str) +
            `</code></pre>`
          );
        },
      });

      // 使用插件
      md.use(MarkdownItKatex);
      md.use(MarkdownItToc, {
        containerClass: 'toc',
        placeholder: '[[toc]]',
      });

      // 处理目录占位符
      let content = this.content.replace(/\[\[toc\]\]/, '[toc]');

      // 渲染为 HTML
      const fullHtml = md.render(content);

      // 提取目录 HTML（支持 <nav class="toc"> 或 <ul class="toc">）
      const tocMatch =
        fullHtml.match(/<nav[^>]*class="toc"[^>]*>[\s\S]*?<\/nav>/i) ||
        fullHtml.match(/<ul[^>]*class="toc"[^>]*>[\s\S]*?<\/ul>/i);

      this.tocHtml = tocMatch ? tocMatch[0] : '';

      // 移除目录部分，保留正文
      this.html = tocMatch ? fullHtml.replace(tocMatch[0], '') : fullHtml;
    },
  },
};
</script>

<!-- 全局样式：作用于 v-html 渲染的内容 -->
<style>
/* 表格样式（必须全局） */
.markdown-body table {
  width: 100%;
  border-collapse: collapse;
  margin: 16px 0;
  font-size: 1em;
}

.markdown-body th,
.markdown-body td {
  padding: 8px 12px;
  border: 1px solid #ddd;
  text-align: left;
  vertical-align: top;
}

.markdown-body th {
  background-color: #f5f7fa;
  font-weight: 600;
}

.markdown-body tr:nth-child(even) {
  background-color: #fafafa;
}

/* KaTeX 公式样式由外部 katex.min.css 提供 */
/* 请确保在 main.js 中引入：import 'katex/dist/katex.min.css' */

/* highlight.js 样式由外部提供 */
/* 请确保在 main.js 中引入：import 'highlight.js/styles/github.css' */
</style>

<!-- 作用于组件自身结构的样式 -->
<style scoped>
.markdown-viewer {
  margin-top: 20px;
}

.markdown-body {
  font-size: 16px;
  color: #333;
  line-height: 1.6;
  max-width: 900px;
  margin: 0 auto;
  padding: 20px;
  background: #fff;
  border: 1px solid #ddd;
  border-radius: 8px;
  box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
}

.markdown-body img {
  max-width: 100%;
  height: auto;
  border-radius: 6px;
  margin: 10px 0;
}

.toc {
  background-color: #f8f9fa;
  border: 1px solid #ddd;
  border-radius: 6px;
  padding: 12px;
}
</style>
