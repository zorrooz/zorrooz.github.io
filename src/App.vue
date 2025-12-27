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
  min-height: calc(100vh - 80px - 80px);
}

.main-content>* {
  flex: 1;
  display: flex;
  flex-direction: column;
}

.main-content .container,
.main-content .container-fluid {
  flex: 1;
  display: flex;
  flex-direction: column;
}

.main-content>*>* {
  flex: 1;
}

.main-content .view-container {
  flex: 1;
  display: flex;
  flex-direction: column;
}

.main-content .view-container>* {
  flex: 1;
}
</style>
