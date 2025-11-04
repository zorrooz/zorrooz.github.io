<!-- AppHeader.vue -->
<template>
  <!-- 外层负责定位和全宽 -->
  <header class="shadow-sm" style="background-color: var(--app-header-bg);">
    <div class="container px-0">
      <nav class="navbar navbar-expand-lg">
        <div class="container-fluid d-flex">
          <!-- Logo -->
          <RouterLink class="navbar-brand" to="/" @mouseover="enlargeLogo" @mouseleave="resetLogo" :style="logoStyle"
            @click="mobileMenuOpen = false">
            <img src="@/assets/icons/gblog.svg" alt="gblog" class="logo-icon" height="40">
          </RouterLink>

          <!-- ========== 移动端按钮 ========== -->
          <div class="d-flex d-lg-none ms-auto">
            <button class="btn btn-sm mx-1 p-2 btn-icon" @click="toggleTheme" @focus="$event.target.blur()">
              <img src="@/assets/icons/change-theme.png" alt="主题" width="20" height="20">
            </button>
            <button class="btn btn-sm mx-1 p-2 btn-icon" @click="toggleLanguage" @focus="$event.target.blur()">
              <img src="@/assets/icons/change-language.png" alt="切换语言" width="20" height="20">
            </button>
            <button class="btn btn-sm mx-1 p-2 btn-icon" @click="onMobileMenuClick"
              style="font-size: 1.2em; background: none; border: none;" @focus="$event.target.blur()">
              <img src="@/assets/icons/menu.png" alt="菜单" width="20" height="20">
            </button>
          </div>

          <!-- ========== 导航菜单（移动端折叠内容）========== -->
          <div :class="['navbar-collapse collapse', { 'show': mobileMenuOpen }]">
            <ul class="navbar-nav mb-2 mb-lg-0">
              <li class="nav-item" v-for="item in navItems" :key="item.text">
                <RouterLink class="nav-link px-3 py-2 rounded d-flex align-items-center" :to="item.href"
                  @click="mobileMenuOpen = false" :class="{ 'active': $route.path === item.href }">
                  <i :class="['fas', item.icon]" style="flex-shrink: 0; width: 20px; margin-right: 8px;"></i>
                  {{ item.text }}
                </RouterLink>
              </li>
            </ul>
          </div>

          <!-- ========== 桌面端按钮 ========== -->
          <div class="d-none d-lg-flex ms-auto">
            <button class="btn btn-sm me-2 btn-icon" @click="toggleTheme" @focus="$event.target.blur()">
              <img src="@/assets/icons/change-theme.png" alt="主题" width="20" height="20">
            </button>
            <button class="btn btn-sm btn-icon" @click="toggleLanguage" @focus="$event.target.blur()">
              <img src="@/assets/icons/change-language.png" alt="切换语言" width="20" height="20">
            </button>
          </div>
        </div>
      </nav>
    </div>
  </header>

  <!-- 移动端统一左侧抽屉：全站风格一致；文章页额外显示目录 -->
  <div v-if="showMobileSidebar" class="mobile-offcanvas d-lg-none" @click.self="closeMobileSidebar">
    <div class="offcanvas-panel border-end rounded-0" :style="{ 'border-color': 'var(--app-border)', 'box-shadow': 'var(--app-offcanvas-shadow)' }">
      <!-- <div class="offcanvas-header d-flex align-items-center justify-content-between">
        <div></div>
        <button class="btn btn-sm p-2 btn-icon" @click="closeMobileSidebar" @focus="$event.target.blur()" aria-label="关闭">
          <i class="bi bi-x-lg"></i>
        </button>
      </div> -->

      <div class="offcanvas-section">
        <div class="offcanvas-card">
          <ul class="list-unstyled m-0">
            <li v-for="item in navItems" :key="item.text" class="my-1">
              <RouterLink
                :to="item.href"
                class="offcanvas-link d-flex align-items-center"
                :class="{ active: $route.path === item.href }"
                @click="closeMobileSidebar"
              >
                <i :class="['fas', item.icon]" style="flex-shrink: 0; width: 20px; margin-right: 8px;"></i>
                <span>{{ item.text }}</span>
              </RouterLink>
            </li>
          </ul>
        </div>
      </div>

      <div v-if="isArticle" class="offcanvas-section">
        <div class="offcanvas-tree offcanvas-card" @click="handleDirectoryClick">
          <NavigationTree />
        </div>
      </div>

    </div>
    <div class="offcanvas-backdrop" @click="closeMobileSidebar"></div>
  </div>
</template>

<script setup>
import { ref, onMounted, onBeforeUnmount, computed, watch } from 'vue';
import { useRoute } from 'vue-router';
import { useI18n } from 'vue-i18n';
import { useAppStore } from '@/stores/app';
import NavigationTree from '@/components/layout/NavigationTree.vue';

const route = useRoute();
const { t, locale } = useI18n();
const appStore = useAppStore();
const mobileMenuOpen = ref(false);
const showMobileSidebar = ref(false);
const logoStyle = ref({});

const navItems = ref([
  { icon: 'fa-layer-group', text: t('categories'), href: '/category' },
  { icon: 'fa-folder-open', text: t('resources'), href: '/resource' },
  { icon: 'fa-info-circle', text: t('about'), href: '/about' },
]);

// 监听语言变化，更新导航项文本
watch(locale, () => {
  navItems.value = [
    { icon: 'fa-layer-group', text: t('categories'), href: '/category' },
    { icon: 'fa-folder-open', text: t('resources'), href: '/resource' },
    { icon: 'fa-info-circle', text: t('about'), href: '/about' },
  ];
});

const isArticle = computed(() => route.name === 'Article');

const enlargeLogo = () => {
  logoStyle.value = {
    transform: 'scale(1.1)',
    transition: 'transform 0.2s ease'
  };
};

const resetLogo = () => {
  logoStyle.value = {};
};

const toggleTheme = () => {
  appStore.toggleTheme();
};

const toggleLanguage = () => {
  appStore.toggleLanguage();
};

const onMobileMenuClick = () => {
  const isMobile = window.innerWidth < 992;
  if (isMobile) {
    openMobileSidebar();
  } else {
    mobileMenuOpen.value = !mobileMenuOpen.value;
  }
};

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

const handleDirectoryClick = (event) => {
  if (event.target.closest('a')) {
    closeMobileSidebar();
  }
};

let _openSidebarHandler = null;
onMounted(() => {
  appStore.initTheme();
  appStore.initLocale();

  _openSidebarHandler = () => openMobileSidebar();
  window.addEventListener('open-mobile-sidebar', _openSidebarHandler);
});

onBeforeUnmount(() => {
  if (_openSidebarHandler) {
    window.removeEventListener('open-mobile-sidebar', _openSidebarHandler);
    _openSidebarHandler = null;
  }
  unlockScroll();
});
</script>

<style scoped>
:global([data-bs-theme="dark"] .btn-icon img) {
  filter: invert(1);
}

/* Scoped styles for AppHeader */
.header {
  background-color: var(--app-header-bg);
}

/* Styles for nav-link are now global in main.scss */

/* 统一移动端抽屉样式（贴近 Bootstrap） */
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
  top: 0; left: 0; bottom: 0;
  width: min(85vw, 320px);
  background: var(--app-offcanvas-bg);
  box-shadow: var(--app-offcanvas-shadow);
  overflow-y: auto;
  -webkit-overflow-scrolling: touch;
  padding: 0.75rem 0.75rem 1rem;
  z-index: 2;
}
.offcanvas-header {
  margin-bottom: 0.5rem;
  padding-bottom: 0.5rem;
}
.offcanvas-section + .offcanvas-section {
  margin-top: 0.75rem;
}
.section-title {
  font-size: 0.95rem;
  color: var(--app-text-muted);
  margin-bottom: 0.25rem;
  font-weight: 600;
}
.offcanvas-link {
  display: block;
  padding: 0.5rem 0.75rem;
  border-radius: 0.5rem;
  color: var(--app-text-muted);
  text-decoration: none;
  transition: background-color 0.15s ease-in-out, color 0.15s ease-in-out, border-color 0.15s ease-in-out;
  line-height: 1.5;
  background-color: transparent;
  border: none;
}
.offcanvas-link:hover,
.offcanvas-link:focus {
  color: var(--app-primary) !important;
  background-color: var(--app-primary-bg-subtle) !important;
}
.offcanvas-link.active {
  color: var(--app-primary) !important;
  background-color: transparent !important;
  font-weight: 500;
}
.offcanvas-link:focus {
  outline: none;
  box-shadow: none;
}
.offcanvas-link i {
  flex-shrink: 0;
  width: 20px;
  margin-right: 8px;
}
.offcanvas-link:hover i,
.offcanvas-link:focus i,
.offcanvas-link.active i {
  color: var(--app-primary) !important;
}
.offcanvas-tree {
  border-top: 0;
}
.offcanvas-card {
  background-color: var(--app-card-bg);
  border: none;
  border-radius: 0.5rem;
  padding: 0.5rem;
  margin: 1rem 0;
}
</style>