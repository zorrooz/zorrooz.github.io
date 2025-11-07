<template>
  <div class="app-container" style="background-color: var(--app-bg);">
    <AppHeader />

    <main class="main-content" style="background-color: var(--app-bg-tertiary);">
      <router-view />
    </main>

    <AppFooter />
    <BackToTop />
  </div>
</template>

<script>
import AppHeader from '@/components/layout/AppHeader.vue'
import AppFooter from '@/components/layout/AppFooter.vue'
import BackToTop from '@/components/widgets/BackToTop.vue'

export default {
  name: 'App',
  components: {
    AppHeader,
    AppFooter,
    BackToTop
  }
}
</script>

<style scoped>
.app-container {
  display: flex;
  flex-direction: column;
  min-height: 100vh;
}

.main-content {
  flex: 1;
  display: flex;
  flex-direction: column;
  min-height: calc(100vh - 80px - 80px); /* 减去header和footer的预估高度 */
}

/* 确保router-view内容区域能够撑满剩余空间 */
.main-content > * {
  flex: 1;
  display: flex;
  flex-direction: column;
}

/* 针对具体视图的样式调整 */
.main-content .container,
.main-content .container-fluid {
  flex: 1;
  display: flex;
  flex-direction: column;
}

/* 确保所有直接子元素都能撑满 */
.main-content > * > * {
  flex: 1;
}

/* 确保view-container类能够正确撑满 */
.main-content .view-container {
  flex: 1;
  display: flex;
  flex-direction: column;
}

.main-content .view-container > * {
  flex: 1;
}
</style>
