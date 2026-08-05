<!-- NavActions.vue — 顶栏操作按钮组（搜索 / 主题 / 语言 / 移动端菜单），移动与桌面共享一份实现 -->
<script setup lang="ts">
import { computed } from 'vue'
import { useRoute, useRouter } from 'vue-router'
import { useI18n } from 'vue-i18n'
import { useThemeStore } from '@/stores/theme'
import { switchLocale } from '@/utils/localePath'

withDefaults(defineProps<{ mobile?: boolean }>(), { mobile: false })
const emit = defineEmits<{ 'open-search': []; 'toggle-menu': [] }>()

const { t } = useI18n()
const route = useRoute()
const router = useRouter()
const themeStore = useThemeStore()

const isDark = computed(() => {
  if (themeStore.theme === 'dark') return true
  if (themeStore.theme === 'light') return false
  return typeof window !== 'undefined' && window.matchMedia('(prefers-color-scheme: dark)').matches
})

const blurFocus = (event: FocusEvent) => {
  ;(event.currentTarget as HTMLButtonElement | null)?.blur()
}

const toggleTheme = () => themeStore.toggleTheme()
const toggleLanguage = () => switchLocale(router, route)
</script>

<template>
  <div
    :class="[
      'd-flex',
      mobile
        ? 'd-lg-none ms-auto app-nav__actions app-nav__actions--mobile'
        : 'd-none d-lg-flex ms-auto app-nav__actions',
    ]"
  >
    <button
      class="icon-btn"
      @click="emit('open-search')"
      @focus="blurFocus"
      :aria-label="t('search')"
    >
      <svg
        class="app-nav__icon"
        width="18"
        height="18"
        viewBox="0 0 24 24"
        fill="none"
        stroke="currentColor"
        stroke-width="1.75"
        stroke-linecap="round"
        stroke-linejoin="round"
        aria-hidden="true"
      >
        <circle cx="11" cy="11" r="8" />
        <line x1="21" y1="21" x2="16.65" y2="16.65" />
      </svg>
    </button>

    <button class="icon-btn" @click="toggleTheme" @focus="blurFocus" :aria-label="t('theme')">
      <svg
        v-if="isDark"
        class="app-nav__icon"
        width="18"
        height="18"
        viewBox="0 0 24 24"
        fill="none"
        stroke="currentColor"
        stroke-width="1.75"
        stroke-linecap="round"
        stroke-linejoin="round"
        aria-hidden="true"
      >
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
      <svg
        v-else
        class="app-nav__icon"
        width="18"
        height="18"
        viewBox="0 0 24 24"
        fill="none"
        stroke="currentColor"
        stroke-width="1.75"
        stroke-linecap="round"
        stroke-linejoin="round"
        aria-hidden="true"
      >
        <path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z" />
      </svg>
    </button>

    <button class="icon-btn" @click="toggleLanguage" @focus="blurFocus" :aria-label="t('language')">
      <svg
        class="app-nav__icon"
        width="18"
        height="18"
        viewBox="0 0 24 24"
        fill="none"
        stroke="currentColor"
        stroke-width="1.75"
        stroke-linecap="round"
        stroke-linejoin="round"
        aria-hidden="true"
      >
        <path d="m5 8 6 6" />
        <path d="m4 14 6-6 2-3" />
        <path d="M2 5h12" />
        <path d="M7 2h1" />
        <path d="m22 22-5-10-5 10" />
        <path d="M14 18h6" />
      </svg>
    </button>

    <button
      v-if="mobile"
      class="icon-btn"
      @click="emit('toggle-menu')"
      @focus="blurFocus"
      :aria-label="t('menu')"
    >
      <svg
        class="app-nav__icon"
        width="18"
        height="18"
        viewBox="0 0 24 24"
        fill="none"
        stroke="currentColor"
        stroke-width="1.75"
        stroke-linecap="round"
        stroke-linejoin="round"
        aria-hidden="true"
      >
        <line x1="3" y1="12" x2="21" y2="12" />
        <line x1="3" y1="6" x2="21" y2="6" />
        <line x1="3" y1="18" x2="21" y2="18" />
      </svg>
    </button>
  </div>
</template>

<style scoped>
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
</style>
