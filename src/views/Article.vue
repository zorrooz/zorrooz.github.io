<!-- Article.vue -->
<template>
  <div class="container">
    <div class="row justify-content-center px-3">
      <div class="col-12 col-lg-8">
        <!-- <TestFormula /> -->
        <RenderMarkdown :rawMarkdown="rawMarkdown" />
      </div>
    </div>
  </div>
</template>

<script>
import RenderMarkdown from '@/components/layout/RenderMarkdown.vue'
import rawMarkdown from '@/content/demo.md?raw'
// import TestFormula from '@/components/layout/TestFormula.vue'

export default {
  name: 'App',
  components: {
    RenderMarkdown,
    // TestFormula
  },
  data() {
    return {
      rawMarkdown
    }
  }
}
</script>
