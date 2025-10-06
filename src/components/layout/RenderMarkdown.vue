<template>
  <div class="markdown-body" v-html="renderedMarkdown"></div>
</template>

<script>
import { renderMarkdown } from '@/utils/markdownProcessor'

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
  watch: {
    // 监听原始Markdown内容变化，重新渲染
    async rawMarkdown(val) {
      // 1. 等待Markdown解析渲染完成
      this.renderedMarkdown = await renderMarkdown(val)
      // 2. 确保DOM已更新后，发射渲染完成事件
      this.$nextTick(() => {
        this.$emit('markdown-rendered')
      })
    }
  },
  async created() {
    // 初始渲染Markdown
    this.renderedMarkdown = await renderMarkdown(this.rawMarkdown)
    // 初始渲染完成后也发射事件
    this.$nextTick(() => {
      this.$emit('markdown-rendered')
    })
  }
}
</script>

<style>
/* Markdown基础样式 */
.markdown-body {
  box-sizing: border-box;
  width: 100%;
  padding: 0;
}



/* 图片样式保持原样 */
.markdown-body img {
  max-width: 100%;
  height: auto;
  display: block;
  margin: 1em auto;
}

/* 可根据需要添加其他Markdown基础样式 */
.markdown-body p {
  margin-top: 0;
  margin-bottom: 16px;
}

.markdown-body h1,
.markdown-body h2,
.markdown-body h3,
.markdown-body h4,
.markdown-body h5,
.markdown-body h6 {
  margin-top: 24px;
  margin-bottom: 16px;
  font-weight: 600;
  line-height: 1.25;
}
</style>
