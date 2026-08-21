<!-- 顶栏操作按钮组（搜索 / 主题 / 语言 / 移动端菜单），移动与桌面共享一份实现 -->
<script setup lang="ts">
import { computed } from 'vue'
import { useRoute, useRouter } from 'vue-router'
import { useI18n } from 'vue-i18n'
import { useThemeStore } from '@/stores/theme'
import { switchLocale } from '@/utils/navigation'
import IconButton from '@/components/common/IconButton.vue'

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
    <IconButton icon="search" :ariaLabel="t('search')" @click="emit('open-search')" />
    <IconButton
      :icon="isDark ? 'sun' : 'moon'"
      :ariaLabel="t('theme')"
      @click="toggleTheme"
    />
    <IconButton icon="language" :ariaLabel="t('language')" @click="toggleLanguage" />
    <IconButton
      v-if="mobile"
      icon="menu"
      :ariaLabel="t('menu')"
      @click="emit('toggle-menu')"
    />
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
</style>
