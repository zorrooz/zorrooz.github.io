<!-- AppHeader.vue -->
<template>
  <header class="app-header">
    <div class="container px-0">
      <nav class="navbar navbar-expand-lg app-nav">
        <div class="container-fluid d-flex app-nav__inner">
          <RouterLink class="navbar-brand app-nav__brand" :to="toLocalePath('/')" @click="mobileMenuOpen = false">
            <span class="app-nav__wordmark">{{ SITE.author }}<span class="app-nav__apos">’</span>s blog</span>
          </RouterLink>

          <div class="d-flex d-lg-none ms-auto app-nav__actions app-nav__actions--mobile">
            <button class="icon-btn" @click="emit('open-search')" @focus="blurFocus"
              :aria-label="t('search')">
              <svg class="app-nav__icon" width="18" height="18" viewBox="0 0 24 24" fill="none"
                stroke="currentColor" stroke-width="1.75" stroke-linecap="round" stroke-linejoin="round"
                aria-hidden="true">
                <circle cx="11" cy="11" r="8" />
                <line x1="21" y1="21" x2="16.65" y2="16.65" />
              </svg>
            </button>
            <button class="icon-btn" @click="toggleTheme" @focus="blurFocus"
              :aria-label="t('theme')">
              <svg v-if="isDark" class="app-nav__icon" width="18" height="18" viewBox="0 0 24 24"
                fill="none" stroke="currentColor" stroke-width="1.75" stroke-linecap="round"
                stroke-linejoin="round" aria-hidden="true">
                <circle cx="12" cy="12" r="5" />
                <line x1="12" y1="1" x2="12" y2="3" />
                <line x1="12" y1="21" x2="12" y2="23" />
                <line x1="4.22" y1="4.22" x2="5.64" y2="5.64" />
                <line x1="18.36" y1="18.36" x2="19.78" y2="19.78" />
                <line x1="1" y1="12" x2="3" y2="12" />
                <line x1="21" y1="12" x2="23" y2="12" />
                <line x1="4.22" y1="19.78" x2="5.64" y2="18.36" />
                <line x1="18.36" y1="5.64" x2="19.78" y2="4.22" />
              </svg>
              <svg v-else class="app-nav__icon" width="18" height="18" viewBox="0 0 24 24"
                fill="none" stroke="currentColor" stroke-width="1.75" stroke-linecap="round"
                stroke-linejoin="round" aria-hidden="true">
                <path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z" />
              </svg>
            </button>
            <button class="icon-btn" @click="toggleLanguage" @focus="blurFocus"
              :aria-label="t('language')">
              <svg class="app-nav__icon" width="18" height="18" viewBox="0 0 24 24" fill="none"
                stroke="currentColor" stroke-width="1.75" stroke-linecap="round" stroke-linejoin="round"
                aria-hidden="true">
                <path d="m5 8 6 6" />
                <path d="m4 14 6-6 2-3" />
                <path d="M2 5h12" />
                <path d="M7 2h1" />
                <path d="m22 22-5-10-5 10" />
                <path d="M14 18h6" />
              </svg>
            </button>
            <button class="icon-btn" @click="onMobileMenuClick" @focus="blurFocus"
              :aria-label="t('menu')">
              <svg class="app-nav__icon" width="18" height="18" viewBox="0 0 24 24" fill="none"
                stroke="currentColor" stroke-width="1.75" stroke-linecap="round" stroke-linejoin="round"
                aria-hidden="true">
                <line x1="3" y1="12" x2="21" y2="12" />
                <line x1="3" y1="6" x2="21" y2="6" />
                <line x1="3" y1="18" x2="21" y2="18" />
              </svg>
            </button>
          </div>

          <div :class="['navbar-collapse collapse', { 'show': mobileMenuOpen }]">
            <ul class="navbar-nav mb-2 mb-lg-0 app-nav__links">
              <li v-if="navItems.length" class="app-nav__link-divider" aria-hidden="true"></li>
              <li class="nav-item" v-for="item in navItems" :key="item.text">
                <RouterLink class="nav-link app-nav__link" :to="item.href"
                  @click="mobileMenuOpen = false" :class="{ 'active': isActive(item.href) }">
                  {{ item.text }}
                </RouterLink>
              </li>
            </ul>
          </div>

          <div class="d-none d-lg-flex ms-auto app-nav__actions">
            <button class="icon-btn" @click="emit('open-search')" @focus="blurFocus"
              :aria-label="t('search')">
              <svg class="app-nav__icon" width="18" height="18" viewBox="0 0 24 24" fill="none"
                stroke="currentColor" stroke-width="1.75" stroke-linecap="round" stroke-linejoin="round"
                aria-hidden="true">
                <circle cx="11" cy="11" r="8" />
                <line x1="21" y1="21" x2="16.65" y2="16.65" />
              </svg>
            </button>
            <button class="icon-btn" @click="toggleTheme" @focus="blurFocus"
              :aria-label="t('theme')">
              <svg v-if="isDark" class="app-nav__icon" width="18" height="18" viewBox="0 0 24 24"
                fill="none" stroke="currentColor" stroke-width="1.75" stroke-linecap="round"
                stroke-linejoin="round" aria-hidden="true">
                <circle cx="12" cy="12" r="5" />
                <line x1="12" y1="1" x2="12" y2="3" />
                <line x1="12" y1="21" x2="12" y2="23" />
                <line x1="4.22" y1="4.22" x2="5.64" y2="5.64" />
                <line x1="18.36" y1="18.36" x2="19.78" y2="19.78" />
                <line x1="1" y1="12" x2="3" y2="12" />
                <line x1="21" y1="12" x2="23" y2="12" />
                <line x1="4.22" y1="19.78" x2="5.64" y2="18.36" />
                <line x1="18.36" y1="5.64" x2="19.78" y2="4.22" />
              </svg>
              <svg v-else class="app-nav__icon" width="18" height="18" viewBox="0 0 24 24"
                fill="none" stroke="currentColor" stroke-width="1.75" stroke-linecap="round"
                stroke-linejoin="round" aria-hidden="true">
                <path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z" />
              </svg>
            </button>
            <button class="icon-btn" @click="toggleLanguage" @focus="blurFocus"
              :aria-label="t('language')">
              <svg class="app-nav__icon" width="18" height="18" viewBox="0 0 24 24" fill="none"
                stroke="currentColor" stroke-width="1.75" stroke-linecap="round" stroke-linejoin="round"
                aria-hidden="true">
                <path d="m5 8 6 6" />
                <path d="m4 14 6-6 2-3" />
                <path d="M2 5h12" />
                <path d="M7 2h1" />
                <path d="m22 22-5-10-5 10" />
                <path d="M14 18h6" />
              </svg>
            </button>
          </div>
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
              <RouterLink :to="item.href" class="offcanvas-link d-flex align-items-center"
                :class="{ active: isActive(item.href) }" @click="closeMobileSidebar">
                <span>{{ item.text }}</span>
                <i class="fas fa-chevron-right offcanvas-link__chevron"></i>
              </RouterLink>
            </li>
          </ul>
        </div>
      </div>

      <div v-if="isArticle" class="offcanvas-section">
        <div class="offcanvas-head">{{ t('tableOfContents') }}</div>
        <div class="offcanvas-tree offcanvas-card" @click="handleDirectoryClick">
          <NavigationTree />
        </div>
      </div>

    </div>
    <div class="offcanvas-backdrop" @click="closeMobileSidebar"></div>
  </div>
</template>

<script setup lang="ts">
import { ref, onMounted, onBeforeUnmount, computed } from 'vue';
import { useRoute, useRouter } from 'vue-router';
import { useI18n } from 'vue-i18n';
import { useAppStore } from '@/stores/app'
import { SITE } from '@/config/site.ts';
import { toLocalePath } from '@/utils/localePath';
import NavigationTree from '@/components/layout/NavigationTree.vue';

const route = useRoute();
const router = useRouter();
const { t } = useI18n();
const appStore = useAppStore();

const isDark = computed(() => {
  if (appStore.theme === 'dark') return true
  if (appStore.theme === 'light') return false
  return (
    typeof window !== 'undefined' &&
    window.matchMedia('(prefers-color-scheme: dark)').matches
  )
})
const mobileMenuOpen = ref(false);
const showMobileSidebar = ref(false);
const emit = defineEmits(['open-search']);

const blurFocus = (event: FocusEvent) => {
  (event.currentTarget as HTMLButtonElement | null)?.blur()
};

const navItems = computed(() => [
  { text: t('categories'), href: toLocalePath('/category') },
  { text: t('resources'), href: toLocalePath('/resource') },
  { text: t('about'), href: toLocalePath('/about') },
])

const isActive = (href: string) => {
  if (href === route.path) return true
  return href !== toLocalePath('/') && route.path.startsWith(href)
}

const isArticle = computed(() => route.path.includes('/article/'));

const toggleTheme = () => {
  appStore.toggleTheme();
};

const toggleLanguage = () => {
  const m = route.path.match(/^\/(zh|en)/)
  const nextPrefix = m && m[1] === 'en' ? '/zh' : '/en'
  const rest = m ? route.path.slice(m[0].length) : route.path
  router.push({ path: `${nextPrefix}${rest}`, query: route.query })
};

const onMobileMenuClick = () => { const isMobile = window.innerWidth < 992; if (isMobile) openMobileSidebar(); else mobileMenuOpen.value = !mobileMenuOpen.value; };

const lockScroll = () => {
  const docEl = document.documentElement;
  const body = document.body;
  if (docEl) {
    docEl.style.overflow = 'hidden';
    docEl.style.overscrollBehavior = 'contain';
  }
  if (body) {
    body.style.overflow = 'hidden';
    body.style.overscrollBehavior = 'contain';
  }
};

const unlockScroll = () => {
  const docEl = document.documentElement;
  const body = document.body;
  if (docEl) {
    docEl.style.overflow = '';
    docEl.style.overscrollBehavior = '';
  }
  if (body) {
    body.style.overflow = '';
    body.style.overscrollBehavior = '';
  }
};

const openMobileSidebar = () => {
  lockScroll();
  showMobileSidebar.value = true;
};

const closeMobileSidebar = () => {
  showMobileSidebar.value = false;
  unlockScroll();
};

const handleDirectoryClick = (event: MouseEvent) => {
  const target = event.target as Element | null
  if (target && target.closest('a')) {
    closeMobileSidebar();
  }
};

onMounted(() => {
  appStore.initTheme();
  appStore.initLocale();
});

onBeforeUnmount(() => {
  unlockScroll();
});
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

.app-nav__actions {
  display: flex;
  gap: 4px;
  align-items: center;
}

.app-nav__actions--mobile {
  gap: 2px;
}

.app-nav__icon {
  display: block;
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
  transition: background-color 0.14s ease, color 0.14s ease;
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
