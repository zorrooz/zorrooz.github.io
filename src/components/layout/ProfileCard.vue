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

<script>
/* 
  ProfileCard 
  - 显示作者和统计信息，基于 posts/tags 数据
*/
import { useI18n } from 'vue-i18n'
import { loadPosts, loadTags } from '@/utils/contentLoader'

export default {
  name: 'ProfileCard',
  setup() {
    const { t, locale } = useI18n()
    return { t, locale }
  },
  data() {
    return {
      posts: [],
      tags: []
    }
  },
  async created() {
    await this.loadData()
  },
  watch: {
    locale() {
      this.loadData()
    }
  },
  methods: {
    async loadData() {
      try {
        const [postsData, tagsData] = await Promise.all([
          loadPosts(),
          loadTags()
        ]);
        this.posts = postsData || [];
        this.tags = tagsData || [];
      } catch (error) {
        console.error('Failed to load data:', error);
        this.posts = [];
        this.tags = [];
      }
    }
  },
  computed: {
    postCount() {
      return Array.isArray(this.posts) ? this.posts.length : 0
    },
    tagCount() {
      return Array.isArray(this.tags) ? this.tags.length : 0
    },
    totalWords() {
      if (!Array.isArray(this.posts)) return 0
      return this.posts.reduce((sum, p) => {
        const n = typeof p?.wordCount === 'number' ? p.wordCount : 0
        return sum + (Number.isFinite(n) ? n : 0)
      }, 0)
    },
    totalWordsDisplay() {
      const n = this.totalWords
      if (n >= 1_000_000) return (n / 1_000_000).toFixed(n % 1_000_000 ? 1 : 0) + 'M'
      if (n >= 1_000) return (n / 1_000).toFixed(n % 1_000 ? 1 : 0) + 'K'
      return String(n)
    }
  }
}
</script>

<style scoped></style>
