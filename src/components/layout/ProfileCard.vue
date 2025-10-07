<!-- ProfileCard.vue -->
<template>
  <div class="card shadow-sm border-0 bg-white mb-3">
    <!-- 头像区域 -->
    <div class="text-center px-4 pt-4 pb-0">
      <div class="bg-light rounded-circle d-inline-flex align-items-center justify-content-center"
        style="width: 80px; height: 80px;">
        <span class="text-muted">头像</span>
      </div>
    </div>

    <!-- 个人信息区域 -->
    <div class="card-body p-4 text-center typography-body">
      <h3 class="card-title mb-1 fw-bold">zorrooz</h3>
      <p class="card-text text-muted mb-4">
        开发者
      </p>

      <!-- 统计数据 -->
      <div class="row g-0 text-center">
        <div class="col border-end">
          <div class="fw-bold">{{ postCount }}</div>
          <div class="text-muted">文章</div>
        </div>
        <div class="col border-end">
          <div class="fw-bold">{{ tagCount }}</div>
          <div class="text-muted">标签</div>
        </div>
        <div class="col">
          <div class="fw-bold">{{ totalWordsDisplay }}</div>
          <div class="text-muted">字数</div>
        </div>
      </div>
    </div>
  </div>
</template>

<script>
import posts from '@/content/posts.json'
import tags from '@/content/tags.json'

export default {
  name: 'ProfileCard',
  data() {
    return { posts, tags }
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
