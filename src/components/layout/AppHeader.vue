<!-- AppHeader.vue -->
<template>
  <!-- 外层负责定位和全宽 -->
  <header class="site-header shadow-sm">
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
            <button class="btn btn-sm mx-1 p-2 btn-icon" @click="mobileMenuOpen = !mobileMenuOpen"
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
</template>

<script>
export default {
  name: 'AppHeader',
  data() {
    return {
      mobileMenuOpen: false,
      navItems: [
        { icon: 'fa-layer-group', text: '分类', href: '/category' },
        { icon: 'fa-folder-open', text: '资源', href: '/resource' },
        { icon: 'fa-info-circle', text: '关于', href: '/about' },
      ],
      logoStyle: {}
    }
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
      const html = document.documentElement
      const current = html.getAttribute('data-bs-theme') || 'light'
      const next = current === 'light' ? 'dark' : 'light'
      html.setAttribute('data-bs-theme', next)
      localStorage.setItem('theme', next)
    },
    toggleLanguage() {
      const current = this.$i18n.locale.value || 'zh-CN'
      const next = current === 'zh-CN' ? 'en-US' : 'zh-CN'
      this.$i18n.locale.value = next
      localStorage.setItem('locale', next)
      document.documentElement.lang = next
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
  background-color: transparent !important; /* 交由主题文件控制 */
}

.btn-icon:hover,
.btn-icon:focus {
  background-color: transparent !important;
  color: inherit !important;
}
</style>
