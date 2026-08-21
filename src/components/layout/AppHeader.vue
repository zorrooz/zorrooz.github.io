<template>
  <header class="app-header">
    <div class="container px-0">
      <nav class="navbar navbar-expand-lg app-nav">
        <div class="container-fluid d-flex app-nav__inner">
          <RouterLink
            class="navbar-brand app-nav__brand"
            :to="toLocalePath('/')"
            @click="mobileMenuOpen = false"
          >
            <span class="app-nav__wordmark"
              >{{ SITE.author }}<span class="app-nav__apos">’</span>s blog</span
            >
          </RouterLink>

          <NavActions mobile @open-search="emit('open-search')" @toggle-menu="onMobileMenuClick" />

          <div :class="['navbar-collapse collapse', { show: mobileMenuOpen }]">
            <ul class="navbar-nav mb-2 mb-lg-0 app-nav__links">
              <li v-if="navItems.length" class="app-nav__link-divider" aria-hidden="true"></li>
              <li class="nav-item" v-for="item in navItems" :key="item.text">
                <RouterLink
                  class="nav-link app-nav__link"
                  :to="item.href"
                  @click="mobileMenuOpen = false"
                  :class="{ active: isActive(item.href) }"
                >
                  {{ item.text }}
                </RouterLink>
              </li>
            </ul>
          </div>

          <NavActions @open-search="emit('open-search')" />
        </div>
      </nav>
    </div>
  </header>

  <div v-if="showMobileSidebar" class="mobile-offcanvas d-lg-none" @click.self="closeMobileSidebar">
    <div class="offcanvas-panel">
      <div class="offcanvas-section">
        <div class="offcanvas-head">{{ t('menu') }}</div>
        <div class="offcanvas-card">
          <ul class="list-unstyled m-0">
            <li v-for="item in navItems" :key="item.text" class="my-1">
              <RouterLink
                :to="item.href"
                class="offcanvas-link d-flex align-items-center"
                :class="{ active: isActive(item.href) }"
                @click="closeMobileSidebar"
              >
                <span>{{ item.text }}</span>
                <i class="fas fa-chevron-right offcanvas-link__chevron"></i>
              </RouterLink>
            </li>
          </ul>
        </div>
      </div>

      <div v-if="isArticle" class="offcanvas-section">
        <div class="offcanvas-head">{{ t('articleNav') }}</div>
        <div class="offcanvas-tree offcanvas-card" @click="handleDirectoryClick">
          <NavigationTree />
        </div>
      </div>
    </div>
    <div class="offcanvas-backdrop" @click="closeMobileSidebar"></div>
  </div>
</template>

<script setup lang="ts">
import { ref, onMounted, onBeforeUnmount, computed } from 'vue'
import { useRoute } from 'vue-router'
import { useI18n } from 'vue-i18n'
import { useThemeStore } from '@/stores/theme'
import { useLocaleStore } from '@/stores/locale'
import { ARTICLE_ROUTE_PREFIX, SITE } from '@/config'
import { toLocalePath } from '@/utils/navigation'
import { lockScrollOverflow, unlockScrollOverflow } from '@/utils/scroll'
import NavigationTree from '@/components/layout/NavigationTree.vue'
import NavActions from '@/components/layout/NavActions.vue'

const route = useRoute()
const { t } = useI18n()
const themeStore = useThemeStore()
const localeStore = useLocaleStore()

const mobileMenuOpen = ref(false)
const showMobileSidebar = ref(false)
const emit = defineEmits(['open-search'])

const navItems = computed(() => [
  { text: t('categories'), href: toLocalePath('/category') },
  { text: t('resources'), href: toLocalePath('/resource') },
  { text: t('about'), href: toLocalePath('/about') },
])

const isActive = (href: string) => {
  if (href === route.path) return true
  return href !== toLocalePath('/') && route.path.startsWith(href)
}

const isArticle = computed(() => route.path.includes(`${ARTICLE_ROUTE_PREFIX}/`))

const onMobileMenuClick = () => {
  const isMobile = window.innerWidth < 992
  if (isMobile) openMobileSidebar()
  else mobileMenuOpen.value = !mobileMenuOpen.value
}

const openMobileSidebar = () => {
  lockScrollOverflow()
  showMobileSidebar.value = true
}

const closeMobileSidebar = () => {
  showMobileSidebar.value = false
  unlockScrollOverflow()
}

const handleDirectoryClick = (event: MouseEvent) => {
  const target = event.target as Element | null
  if (target && target.closest('a')) {
    closeMobileSidebar()
  }
}

onMounted(() => {
  themeStore.initTheme()
  localeStore.initLocale()
})

onBeforeUnmount(() => {
  unlockScrollOverflow()
})
</script>

<style scoped>
.app-header {
  position: sticky;
  top: 0;
  z-index: 1020;
  background: var(--app-header-bg);
  backdrop-filter: blur(14px) saturate(1.4);
  -webkit-backdrop-filter: blur(14px) saturate(1.4);
  border-bottom: 1px solid var(--line);
}

/* 175%+ 显示缩放下等效视口落入 <1200px：header 内容贴边，不缩在 Bootstrap container 居中 */
@media (max-width: 1199.98px) {
  .app-header .container {
    max-width: 100%;
  }
}

.app-nav {
  height: 64px;
  padding: 0;
}

.app-nav__inner {
  gap: var(--sp-4);
  height: 100%;
}

@media (max-width: 767px) {
  .app-nav__inner {
    padding-inline: 20px;
  }
}

.app-nav__brand {
  display: flex;
  align-items: center;
  gap: 10px;
  font-size: 18px;
  font-weight: 700;
  letter-spacing: -0.01em;
  font-family: var(--font-serif);
  color: var(--fg);
  padding: 0;
  white-space: nowrap;
}

.app-nav__wordmark {
  color: var(--fg);
  line-height: 1;
  padding-top: 1px;
}

/* 弯撇必须用西文字体渲染：思源宋体的 U+2019 是全角引号（1000/1000em），会撑开间距 */
.app-nav__apos {
  font-family: Georgia, 'Times New Roman', 'SourceHanSerifCN', serif;
}

.app-nav__links {
  display: flex;
  height: 64px;
  align-items: center;
  gap: var(--sp-3);
  margin-left: var(--sp-5);
}

.app-nav__links .nav-item {
  display: flex;
  align-items: center;
}

.app-nav__link-divider {
  width: 1px;
  height: 20px;
  background: var(--line-strong);
  margin-right: var(--sp-3);
  flex-shrink: 0;
}

.app-nav__link {
  position: relative;
  display: flex;
  align-items: center;
  height: 100%;
  padding: 0 4px;
  font-size: var(--text-base);
  font-weight: 500;
  color: var(--fg-2);
  transition: color var(--dur-fast) ease;
  letter-spacing: 0.01em;
}

.app-nav__link:hover {
  color: var(--fg);
}

.app-nav__link.active,
.app-nav__link.active:hover,
.app-nav__link.router-link-active,
.app-nav__link.router-link-active:hover,
.navbar .app-nav__link.active,
.navbar .app-nav__link.active:hover,
.navbar .app-nav__link.router-link-active,
.navbar .app-nav__link.router-link-active:hover {
  color: var(--primary);
  font-weight: 700;
}

.mobile-offcanvas {
  position: fixed;
  inset: 0;
  z-index: 1050;
}

.offcanvas-backdrop {
  position: absolute;
  inset: 0;
  background: var(--app-backdrop-bg);
  z-index: 1;
}

.offcanvas-panel {
  position: absolute;
  top: 0;
  left: 0;
  bottom: 0;
  width: min(84vw, 320px);
  background: var(--surface);
  border-right: 1px solid var(--line);
  box-shadow: var(--shadow-lift);
  overflow-y: auto;
  -webkit-overflow-scrolling: touch;
  padding: 1.25rem 1rem 1.5rem;
  z-index: 2;
}

.offcanvas-head {
  font-size: 12px;
  font-weight: 600;
  letter-spacing: 0.08em;
  text-transform: uppercase;
  color: var(--fg-3);
  padding: 0 var(--sp-2) var(--sp-2);
}

.offcanvas-link__chevron {
  margin-left: auto;
  font-size: 10px;
  opacity: 0;
  transition: opacity 0.14s ease;
}

.offcanvas-link:hover .offcanvas-link__chevron,
.offcanvas-link.active .offcanvas-link__chevron {
  opacity: 1;
}

.offcanvas-section + .offcanvas-section {
  margin-top: 0.75rem;
}

.offcanvas-link {
  display: flex;
  align-items: center;
  padding: 0.625rem 0.75rem;
  border-radius: var(--radius-sm);
  color: var(--fg-2);
  text-decoration: none;
  font-size: 14px;
  font-weight: 500;
  transition:
    background-color 0.14s ease,
    color 0.14s ease;
  line-height: 1.5;
  background-color: transparent;
  border: none;
}

.offcanvas-link:hover,
.offcanvas-link:focus {
  color: var(--primary);
  background-color: var(--tint);
}

.offcanvas-link.active {
  color: var(--primary);
  background-color: var(--tint);
  font-weight: 600;
}

.offcanvas-link:focus {
  outline: none;
  box-shadow: none;
}

.offcanvas-link i {
  flex-shrink: 0;
  font-size: 12px;
}

.offcanvas-link:hover i,
.offcanvas-link:focus i,
.offcanvas-link.active i {
  color: var(--primary);
}

.offcanvas-tree {
  border-top: 0;
}

.offcanvas-card {
  background-color: var(--surface);
  border: 1px solid var(--line);
  border-radius: var(--radius);
  padding: 0.5rem;
  margin: 0.5rem 0;
}
</style>
