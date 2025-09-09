<!-- Article.vue -->
<template>
  <div class="container">
    <div class="row justify-content-center">
      <div class="col-12 col-lg-8">
        <RenderMarkdown :rawMarkdown="rawMarkdown" />
      </div>
    </div>
  </div>
</template>

<script>
import RenderMarkdown from '@/components/layout/RenderMarkdown.vue'
import rawMarkdown from '@/content/demo.md?raw'

export default {
  name: 'App',
  components: {
    RenderMarkdown
  },
  data() {
    return {
      rawMarkdown
    }
  }
}
</script>
