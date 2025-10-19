<!-- RecentPosts.vue -->
<template>
  <div class="card shadow-sm border-0" style="max-width: 280px; background-color: var(--app-card-bg);">
    <div class="card-header border-0 pt-3 px-3 pb-2" style="background-color: var(--app-card-bg);">
      <h6 class="m-0 small" style="color: var(--app-text-muted);">最近文章</h6>
    </div>
    <ul class="list-group list-group-flush">
      <li class="list-group-item px-3 py-2" v-for="post in recentPosts" :key="post.id">
        <div class="d-flex flex-column">
          <small style="color: var(--app-text-muted); margin-bottom: 0.25rem;">{{ post.date }}</small>
          <span style="color: var(--app-text);">{{ post.title }}</span>
        </div>
      </li>
    </ul>
  </div>
</template>

<script>
export default {
  name: 'RecentPosts',
  data() {
    return {
      recentPosts: [
        { id: 1, title: 'Vue 3 Composition API详解', date: '今天' },
        // ...其他文章
      ]
    }
  },
  props: {
    posts: {
      type: Array,
      default: () => []
    }
  },
  mounted() {
    if (this.posts.length > 0) {
      this.recentPosts = this.posts;
    }
  }
}
</script>

<style scoped>
/* 无自定义CSS */
</style>
