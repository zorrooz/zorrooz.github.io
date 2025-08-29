<template>
  <div class="language-switcher">
    <button
      @click="setLocale('zh-CN')"
      :class="{ active: current === 'zh-CN' }"
    >
      中文
    </button>
    <button
      @click="setLocale('en-US')"
      :class="{ active: current === 'en-US' }"
    >
      English
    </button>
  </div>
</template>

<script setup>
import { ref } from 'vue'
import { useI18n } from 'vue-i18n'

const { locale } = useI18n()
const current = ref(locale.value)

// 切换语言
const setLocale = (lang) => {
  locale.value = lang
  current.value = lang
  localStorage.setItem('locale', lang)
  document.documentElement.lang = lang
}
</script>

<style scoped>
.language-switcher button {
  margin: 5px;
  padding: 8px 12px;
  border: 1px solid #ccc;
  background: #f8f8f8;
  cursor: pointer;
}

.language-switcher button.active {
  background: #0d6efd;
  color: white;
  border-color: #0d6efd;
}
</style>
