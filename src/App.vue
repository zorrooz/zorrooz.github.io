<template>
  <div class="app-container bg-light">
    <AppHeader />

    <main class="main-content">
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
<style scoped></style>
