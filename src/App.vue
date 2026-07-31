<template>
  <div class="app-container" style="background-color: var(--app-bg);">
    <AppHeader @open-search="searchOpen = true" />

    <main class="main-content" style="background-color: var(--app-bg-tertiary);">
      <router-view />
    </main>

    <AppFooter />
    <BackToTop />
    <SearchModal v-if="searchOpen" @close="searchOpen = false" />
  </div>
</template>

<script setup lang="ts">
import { defineAsyncComponent, onMounted, ref, watch } from 'vue'
import { useRoute } from 'vue-router'
import { useAppStore } from '@/stores/app'
import AppHeader from '@/components/layout/AppHeader.vue'
import AppFooter from '@/components/layout/AppFooter.vue'
import BackToTop from '@/components/widgets/BackToTop.vue'

const searchOpen = ref(false)
const SearchModal = defineAsyncComponent(() => import('@/components/widgets/SearchModal.vue'))

const route = useRoute()
const appStore = useAppStore()

const syncLocaleFromPath = (path: string) => {
  const m = path.match(/^\/(zh|en)(\/|$)/)
  if (m) appStore.setLocale(m[1] === 'en' ? 'en-US' : 'zh-CN')
}

watch(() => route.fullPath, syncLocaleFromPath)
onMounted(() => syncLocaleFromPath(route.fullPath))
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
