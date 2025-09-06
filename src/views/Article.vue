<!-- Article.vue -->
<template>
  <div class="container py-5">
    <div class="row justify-content-center">
      <div class="col-lg-2 d-none d-lg-block">
        <NavigationTree />
      </div>
      <div class="col-lg-10 col-xl-8">
        你好啊
        <Markdown :content="markdownContent" />
      </div>
      <div class="col-lg-2 d-none d-lg-block">
        <OnThisPage />
      </div>
    </div>
  </div>
</template>

<script>
import NavigationTree from '@/components/layout/NavigationTree.vue'
import OnThisPage from '@/components/layout/OnThisPage.vue'
import Markdown from '@/components/layout/Markdown.vue'
import markdownContent from '@/content/demo.md?raw'

export default {
  name: 'ArticleView',
  components: {
    NavigationTree,
    OnThisPage,
    Markdown,
  },
  data() {
    return {
      markdownContent
    }
  }
}
</script>

<style scoped></style>
