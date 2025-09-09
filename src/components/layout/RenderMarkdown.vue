<!-- RenderMarkdown.vue -->
<template>
  <div class="markdown-body py-4">
    <div v-html="renderedMarkdown"></div>
  </div>
</template>

<script>
import MarkdownIt from 'markdown-it'

export default {
  props: ['rawMarkdown'],
  computed: {
    renderedMarkdown() {
      const md = new MarkdownIt()
      return md.render(this.rawMarkdown)
    }
  }
}
</script>

<style scoped></style>
