<template>
  <div class="page-wrap">
    <AppHeader @open-search="searchOpen = true" />

    <main class="main-content">
      <router-view v-slot="{ Component }">
        <transition name="view" mode="out-in">
          <component :is="Component" />
        </transition>
      </router-view>
    </main>

    <AppFooter />
    <BackToTop />
    <SearchModal v-if="searchOpen" @close="searchOpen = false" />
  </div>
</template>

<script setup lang="ts">
import { computed, defineAsyncComponent, onBeforeUnmount, onMounted, ref, watch } from 'vue'
import { useRoute } from 'vue-router'
import { useI18n } from 'vue-i18n'
import { useHead } from '@unhead/vue'
import { useAppStore } from '@/stores/app'
import { localeFromPrefix, SITE } from '@/config/site.ts'
import AppHeader from '@/components/layout/AppHeader.vue'
import AppFooter from '@/components/layout/AppFooter.vue'
import BackToTop from '@/components/widgets/BackToTop.vue'

const searchOpen = ref(false)
const SearchModal = defineAsyncComponent(() => import('@/components/widgets/SearchModal.vue'))

const route = useRoute()
const appStore = useAppStore()
const { locale } = useI18n()

// html lang 与当前语言同步（SSG 产物与运行时一致，保证字体/拼写按语言环境渲染）
// SSG 渲染时 route.path 已含 /zh、/en 前缀，运行时 hash 路由不含——统一剥离后拼接
const pagePath = computed(() => {
  const p = route.path.replace(/^\/(zh|en)(?=\/|$)/, '')
  return p === '/' || p === '' ? '' : p
})
const localePrefix = computed(() => (locale.value === 'en-US' ? '/en' : '/zh'))

useHead({
  htmlAttrs: {
    lang: () => (locale.value === 'en-US' ? 'en' : 'zh-CN'),
  },
  // SEO：canonical + 双语 hreflang（SSG 产物为 /zh、/en 物理路径，避免重复内容）
  // 文章页 en 路径带 -en 后缀，alternate 需做后缀映射
  link: computed(() => {
    const isArticle = pagePath.value.includes('/article/')
    const zhAlt = isArticle ? pagePath.value.replace(/-en$/, '') : pagePath.value
    const enAlt = isArticle && !pagePath.value.endsWith('-en') ? `${pagePath.value}-en` : pagePath.value
    return [
      { rel: 'canonical', href: `${SITE.url}${localePrefix.value}${pagePath.value}` },
      { rel: 'alternate', hreflang: 'zh-CN', href: `${SITE.url}/zh${zhAlt}` },
      { rel: 'alternate', hreflang: 'en', href: `${SITE.url}/en${enAlt}` },
      { rel: 'alternate', hreflang: 'x-default', href: `${SITE.url}/zh${zhAlt}` },
    ]
  }),
  meta: computed(() => [
    { property: 'og:url', content: `${SITE.url}${localePrefix.value}${pagePath.value}` },
  ]),
})

const syncLocaleFromPath = (path: string) => {
  const m = path.match(/^\/(zh|en)(\/|$)/)
  if (m) appStore.setLocale(localeFromPrefix(m[1]))
}

const isTypingTarget = (target: EventTarget | null) => {
  if (!(target instanceof HTMLElement)) return false
  return target.tagName === 'INPUT' || target.tagName === 'TEXTAREA' || target.isContentEditable
}

const onGlobalKeydown = (e: KeyboardEvent) => {
  if (e.key === '/' && !isTypingTarget(e.target)) {
    e.preventDefault()
    searchOpen.value = true
  }
}

watch(() => route.fullPath, syncLocaleFromPath)
onMounted(() => {
  syncLocaleFromPath(route.fullPath)
  window.addEventListener('keydown', onGlobalKeydown)
})
onBeforeUnmount(() => {
  window.removeEventListener('keydown', onGlobalKeydown)
})
</script>

<style scoped>
.main-content {
  flex: 1;
  display: flex;
  flex-direction: column;
  min-height: calc(100vh - 64px - 80px);
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

.view-enter-active {
  transition: opacity 0.22s ease;
}

.view-leave-active {
  transition: opacity 0.14s ease;
}

.view-enter-from {
  opacity: 0;
}

.view-leave-to {
  opacity: 0;
}

@media (prefers-reduced-motion: reduce) {
  .view-enter-active,
  .view-leave-active {
    transition: none;
  }
}
</style>
