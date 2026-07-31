<!-- ProfileCard.vue -->
<template>
  <div class="card shadow-sm border-0 mb-3" :style="{ backgroundColor: 'var(--app-card-bg)' }">
    <div class="text-center px-4 pt-4 pb-0">
      <div class="rounded-circle d-inline-flex align-items-center justify-content-center"
        :style="{ width: '80px', height: '80px', backgroundColor: 'var(--app-bg-light)' }">
        <span :style="{ color: 'var(--app-text-muted)' }">{{ t('avatar') }}</span>
      </div>
    </div>

    <div class="card-body p-4 text-center typography-body">
      <h3 class="card-title mb-1 fw-bold" :style="{ color: 'var(--app-text)' }">zorrooz</h3>
      <p class="card-text mb-4" :style="{ color: 'var(--app-text-muted)' }">
        {{ t('developer') }}
      </p>

      <div class="row g-0 text-center" :style="{ 'border-color': 'var(--app-border)' }">
        <div class="col border-end">
          <div class="fw-bold" :style="{ color: 'var(--app-stat-num-color)' }">{{ postCount }}</div>
          <div :style="{ color: 'var(--app-text-muted)' }">{{ t('articles') }}</div>
        </div>
        <div class="col border-end">
          <div class="fw-bold" :style="{ color: 'var(--app-stat-num-color)' }">{{ tagCount }}</div>
          <div :style="{ color: 'var(--app-text-muted)' }">{{ t('tags') }}</div>
        </div>
        <div class="col">
          <div class="fw-bold" :style="{ color: 'var(--app-stat-num-color)' }">{{ totalWordsDisplay }}</div>
          <div :style="{ color: 'var(--app-text-muted)' }">{{ t('words') }}</div>
        </div>
      </div>
    </div>
  </div>
</template>

<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import { useI18n } from 'vue-i18n'
import { loadPosts, loadTags } from '@/utils/contentLoader'

const { t, locale } = useI18n()

const posts = ref<any[]>([])
const tags = ref<any[]>([])

function loadData() {
  try {
    posts.value = loadPosts() || []
    tags.value = loadTags() || []
  } catch (error) {
    console.error('Failed to load data:', error)
    posts.value = []
    tags.value = []
  }
}

loadData()

watch(locale, () => loadData())

const postCount = computed(() => (Array.isArray(posts.value) ? posts.value.length : 0))
const tagCount = computed(() => (Array.isArray(tags.value) ? tags.value.length : 0))
const totalWords = computed(() => {
  if (!Array.isArray(posts.value)) return 0
  return posts.value.reduce((sum, p) => {
    const n = typeof p?.wordCount === 'number' ? p.wordCount : 0
    return sum + (Number.isFinite(n) ? n : 0)
  }, 0)
})
const totalWordsDisplay = computed(() => {
  const n = totalWords.value
  if (n >= 1_000_000) return (n / 1_000_000).toFixed(n % 1_000_000 ? 1 : 0) + 'M'
  if (n >= 1_000) return (n / 1_000).toFixed(n % 1_000 ? 1 : 0) + 'K'
  return String(n)
})
</script>

<style scoped></style>
