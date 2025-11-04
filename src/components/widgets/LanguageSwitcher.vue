<!-- LanguageSwitcher.vue -->
<template>
  <div class="language-switcher">
    <button class="btn btn-sm btn-icon" @click="toggleLanguage" @focus="$event.target.blur()">
      <img src="@/assets/icons/change-language.png" alt="切换语言" width="20" height="20">
      <span class="ms-1">{{ currentLanguage }}</span>
    </button>
  </div>
</template>

<script setup>
import { computed } from 'vue'
import { useI18n } from 'vue-i18n'
import { useAppStore } from '@/stores/app'

const { t, locale } = useI18n()
const appStore = useAppStore()

const currentLanguage = computed(() => {
  return locale.value === 'zh-CN' ? '中文' : 'English'
})

const toggleLanguage = () => {
  appStore.toggleLanguage()
}
</script>

<style scoped>
.language-switcher {
  display: inline-flex;
  align-items: center;
}

.btn-icon {
  display: inline-flex;
  align-items: center;
  padding: 0.375rem 0.75rem;
  border: 1px solid var(--app-border);
  background: var(--app-card-bg);
  color: var(--app-text);
  transition: all 0.2s ease;
}

.btn-icon:hover {
  border-color: var(--app-primary);
  color: var(--app-primary);
  background: var(--app-primary-bg-subtle);
}

:global([data-bs-theme="dark"] .btn-icon img) {
  filter: invert(1);
}
</style>