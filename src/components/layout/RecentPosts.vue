<!-- RecentPosts.vue -->
<template>
  <div class="card shadow-sm border-0 bg-white" style="max-width: 280px;">
    <div class="card-header bg-white border-0 pt-3 px-3 pb-2">
      <h6 class="text-secondary m-0 small">最近文章</h6>
    </div>
    <ul class="list-group list-group-flush">
      <li class="list-group-item px-3 py-2" v-for="post in recentPosts" :key="post.id">
        <div class="d-flex flex-column">
          <small class="text-muted mb-1">{{ post.date }}</small>
          <span class="text-dark">{{ post.title }}</span>
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
