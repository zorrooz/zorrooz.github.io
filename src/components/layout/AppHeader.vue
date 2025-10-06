<!-- AppHeader.vue -->
<template>
  <!-- 外层负责定位和全宽 -->
  <header class="bg-white shadow-sm">
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
    <div class="offcanvas-panel border-end rounded-0 shadow-sm">
      <div class="offcanvas-header d-flex align-items-center justify-content-between">
        <div></div>
        <button class="btn btn-sm p-2 btn-icon" @click="closeMobileSidebar" @focus="$event.target.blur()" aria-label="关闭">
          <i class="bi bi-x-lg"></i>
        </button>
      </div>

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

<script>
import { setTheme } from '@/utils/themeManager'
import NavigationTree from '@/components/layout/NavigationTree.vue'

export default {
  name: 'AppHeader',
  components: { NavigationTree },
  data() {
    return {
      mobileMenuOpen: false,
      showMobileSidebar: false,
      navItems: [
        { icon: 'fa-layer-group', text: '分类', href: '/category' },
        { icon: 'fa-folder-open', text: '资源', href: '/resource' },
        { icon: 'fa-info-circle', text: '关于', href: '/about' },
      ],
      logoStyle: {},
      theme: 'light'
    }
  },
  computed: {
    isArticle() { return this.$route?.name === 'Article'; }
  },
  methods: {
    enlargeLogo() {
      this.logoStyle = {
        transform: 'scale(1.1)',
        transition: 'transform 0.2s ease'
      };
    },
    resetLogo() {
      this.logoStyle = {};
    },
    toggleTheme() {
      this.theme = this.theme === 'light' ? 'dark' : 'light';
      setTheme(this.theme);
    },
    toggleLanguage() {
      console.log('切换语言');
    },
    onMobileMenuClick() {
      const isMobile = window.innerWidth < 992;
      if (isMobile) {
        this.openMobileSidebar();
      } else {
        this.mobileMenuOpen = !this.mobileMenuOpen;
      }
    },
    openMobileSidebar() {
      this.showMobileSidebar = true;
    },
    closeMobileSidebar() {
      this.showMobileSidebar = false;
    },
    handleDirectoryClick(event) {
      // If a link inside the tree is clicked, close the sidebar.
      if (event.target.closest('a')) {
        this.closeMobileSidebar();
      }
    }
  },
  mounted() {
    const saved = localStorage.getItem('theme');
    if (saved) this.theme = saved;
    // main.js 已在启动时初始化主题，这里仅同步本地状态即可。

    // 监听来自 Article.vue 的移动端抽屉打开事件
    this._openSidebarHandler = () => this.openMobileSidebar();
    window.addEventListener('open-mobile-sidebar', this._openSidebarHandler);
  },
  beforeUnmount() {
    if (this._openSidebarHandler) {
      window.removeEventListener('open-mobile-sidebar', this._openSidebarHandler);
      this._openSidebarHandler = null;
    }
  }
}
</script>

<style scoped>
.nav-link.active {
  color: #047AFF !important;
  background-color: transparent !important;
  font-weight: 500;
}

.nav-link:hover {
  color: #047AFF !important;
  background-color: #E6F0FF !important;
}

.btn-icon:hover,
.btn-icon:focus {
  background-color: transparent !important;
  color: inherit !important;
}
/* 统一移动端抽屉样式（贴近 Bootstrap） */
.mobile-offcanvas {
  position: fixed;
  inset: 0;
  z-index: 1050;
}
.offcanvas-backdrop {
  position: absolute;
  inset: 0;
  background: rgba(0,0,0,0.25);
  z-index: 1;
}
.offcanvas-panel {
  position: absolute;
  top: 0; left: 0; bottom: 0;
  width: min(85vw, 320px);
  background: #e9ecef;
  box-shadow: 0 0.25rem 0.75rem rgba(0,0,0,0.08);
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
  color: var(--bs-gray-700);
  margin-bottom: 0.25rem;
  font-weight: 600;
}
.offcanvas-link {
  display: block;
  padding: 0.5rem 0.75rem;
  border-radius: 0.5rem;
  color: var(--bs-gray-700);
  text-decoration: none;
  transition: background-color 0.15s ease-in-out, color 0.15s ease-in-out, border-color 0.15s ease-in-out;
  line-height: 1.5;
  background-color: transparent;                  /* 由整体卡片承载背景 */
  border: none;
}
.offcanvas-link:hover,
.offcanvas-link:focus {
  color: #047AFF !important;                      /* 与 Header hover 一致 */
  background-color: #E6F0FF !important;           /* 与 Header hover 一致 */
}
.offcanvas-link.active {
  color: #047AFF !important;                      /* 与 Header active 一致 */
  background-color: transparent !important;       /* 与桌面一致：透明背景 */
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
  /* color is inherited */
}
.offcanvas-link:hover i,
.offcanvas-link:focus i,
.offcanvas-link.active i {
  color: #047AFF !important; /* 悬停/激活变蓝，与桌面一致 */
}
.offcanvas-tree {
  border-top: 0;
}
.offcanvas-card {
  background-color: #fff;                         /* 目录卡片：白色圆角矩形 */
  border: none;
  border-radius: 0.5rem;
  padding: 0.5rem;
}
</style>
