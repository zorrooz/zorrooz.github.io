<!-- RenderMarkdown.vue -->
<template>
  <div class="markdown-body py-4" v-html="renderedMarkdown"></div>
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
    async rawMarkdown(val) {
      this.renderedMarkdown = await renderMarkdown(val)
    }
  },
  async created() {
    this.renderedMarkdown = await renderMarkdown(this.rawMarkdown)
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
</style>
