<template>
  <div class="markdown-body py-4" v-html="renderedMarkdown"></div>
</template>

<script>
import { renderMarkdown } from '@/utils/markdownProcessor'
import 'katex/dist/katex.min.css'

export default {
  props: {
    rawMarkdown: {
      type: String,
      default: ''
    }
  },
  data() {
    return {
      renderedMarkdown: ''
    }
  },
  async mounted() {
    this.renderedMarkdown = await renderMarkdown(this.rawMarkdown)
  },
  watch: {
    async rawMarkdown(newVal) {
      this.renderedMarkdown = await renderMarkdown(newVal)
    }
  }
}
</script>

<style>
.markdown-body img {
  max-width: 100%;
  height: auto;
  display: block;
  margin: 1em auto;
}

/* 新增引用块样式 */
.markdown-body blockquote {
  border-left: 4px solid #e2e8f0;
  /* 左侧竖线 */
  padding: 0.5em 1em;
  /* 内边距 */
  margin: 1em 0;
  /* 外边距 */
  color: #4b5563;
  /* 文字颜色 */
  background-color: #f8fafc;
  /* 背景色，可选 */
}

/* 引用块内的段落样式 */
.markdown-body blockquote p {
  margin: 0.5em 0;
  line-height: 1.6;
}
</style>
