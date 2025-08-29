<template>
  <button
    class="btn p-2"
    :class="theme === 'dark' ? 'btn-outline-light' : 'btn-outline-secondary'"
    @click="toggleTheme"
    aria-label="切换主题"
    style="width: 40px; height: 40px; padding: 0; border-radius: 50%; display: flex; align-items: center; justify-content: center;"
  >
    <svg
      v-if="theme === 'light'"
      xmlns="http://www.w3.org/2000/svg"
      width="20"
      height="20"
      fill="currentColor"
      viewBox="0 0 16 16"
    >
      <!-- 暗色模式图标：月亮 -->
      <path d="M6 .278a.768.768 0 0 1 .08.858 7.208 7.208 0 0 0-.878 3.46c0 4.021 3.278 7.277 7.318 7.277.527 0 1.04-.055 1.533-.16a.787.787 0 0 1 .81.686.733.733 0 0 1-.031.293c-.034.214-.121.4.258.4h.059c.442 0 .74.37.742.807A8.05 8.05 0 0 1 13.99 15 8.049 8.049 0 0 1 8 15.989 8.05 8.05 0 0 1 3.011 15a8.05 8.05 0 0 1-.99-14.721z"/>
    </svg>

    <svg
      v-else
      xmlns="http://www.w3.org/2000/svg"
      width="20"
      height="20"
      fill="currentColor"
      viewBox="0 0 16 16"
    >
      <!-- 亮色模式图标：太阳 -->
      <path d="M8 12a4 4 0 1 0 0-8 4 4 0 0 0 0 8zM8 0a.5.5 0 0 1 .5.5v2a.5.5 0 0 1-1 0v-2A.5.5 0 0 1 8 0zm0 13a.5.5 0 0 1 .5.5v2a.5.5 0 0 1-1 0v-2A.5.5 0 0 1 8 13zm7.657-7.657a.5.5 0 0 1 0 .707l-1.414 1.414a.5.5 0 1 1-.707-.707l1.414-1.414a.5.5 0 0 1 .707 0zm-13.214 0a.5.5 0 0 1 .707 0l1.414 1.414a.5.5 0 1 1-.707.707L1.414 7.05a.5.5 0 0 1 0-.707zm9.758 7.071a.5.5 0 0 1 0 .707l-1.414 1.414a.5.5 0 1 1-.707-.707l1.414-1.414a.5.5 0 0 1 .707 0zm-7.071-7.071a.5.5 0 0 1 .707 0l1.414 1.414a.5.5 0 1 1-.707.707L7.05 8.95a.5.5 0 0 1 0-.707zm7.071-7.071a.5.5 0 0 1 .707 0l1.414 1.414a.5.5 0 1 1-.707.707L14.95 4.05a.5.5 0 0 1 0-.707zm-7.071 7.071a.5.5 0 0 1 .707 0l1.414 1.414a.5.5 0 1 1-.707.707l-1.414-1.414a.5.5 0 0 1 0-.707zm-7.071-7.071a.5.5 0 0 1 .707 0l1.414 1.414a.5.5 0 1 1-.707.707L4.05 4.05a.5.5 0 0 1 0-.707z"/>
    </svg>
  </button>
</template>

<script setup>
import { ref, onMounted, watch } from 'vue';

const theme = ref('light');

// 初始化：读取保存的主题
onMounted(() => {
  const saved = localStorage.getItem('theme');
  if (saved) theme.value = saved;
  applyTheme(theme.value);
});

// 切换主题
const toggleTheme = () => {
  theme.value = theme.value === 'light' ? 'dark' : 'light';
};

// 应用主题
const applyTheme = (mode) => {
  document.documentElement.setAttribute('data-bs-theme', mode);
  localStorage.setItem('theme', mode);
};

// 监听变化
watch(theme, applyTheme);
</script>

<style scoped>
/* 可选：添加悬停效果 */
button:hover svg {
  transform: scale(1.1);
  transition: transform 0.2s;
}
</style>
