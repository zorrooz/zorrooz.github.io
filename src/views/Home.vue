<!-- Home.vue -->
<template>
  <div class="container">
    <div class="row py-4 px-0">
      <!-- 中间主内容区域 -->
      <div class="col-12 col-lg-9 order-1 order-lg-2 typography-body" ref="mainContent">
        <div class="row">
          <div class="col">


            <div v-if="currentTag" class="mb-3 d-flex align-items-center gap-2">
              <span>{{ filteredByText }}：</span>
              <span class="current-tag-chip d-inline-flex">
                <span># {{ currentTag }}</span>
                <button class="chip-close" @click="clearTag">×</button>
              </span>
            </div>
            <PostList :docs="filteredDocs" :perPage="5" />
          </div>
        </div>
      </div>

      <!-- 左侧边栏 -->
      <div class="col-12 col-lg-3 order-2 order-lg-1" ref="sidebarContainer">
        <div class="sticky-sidebar" ref="sidebarContent">
          <div class="d-flex flex-column w-100 gap-4">
            <ProfileCard class="w-100" />
            <TagCloud class="w-100" :tagData="tagList" />
          </div>
        </div>
      </div>
    </div>
  </div>
</template>

<script>
import ProfileCard from '@/components/layout/ProfileCard.vue'
import TagCloud from '@/components/layout/TagCloud.vue'
import PostList from '@/components/layout/PostList.vue'
import { useI18n } from 'vue-i18n'
import { loadPosts } from '@/utils/contentLoader'


export default {
  name: 'HomeView',
  setup() {
    const { t, locale } = useI18n()
    return { t, locale }
  },
  components: { ProfileCard, TagCloud, PostList },
  data() {
    return {
      postData: []
    }
  },
  computed: {
    tagsText() {
      return this.t('tags')
    },
    filteredByText() {
      return this.t('filteredBy')
    },
    currentTag() { return this.$route.query.tag || '' },
    filteredDocs() {
      const tag = this.currentTag
      if (!tag) return this.postData
      return this.postData.filter(p => Array.isArray(p.tags) && p.tags.includes(tag))
    },
    tagList() {
      const set = new Set()
      this.postData.forEach(p => (p.tags || []).forEach(t => set.add(t)))
      return Array.from(set).sort()
    }
  },
  async created() {
    await this.loadPostData()
  },
  mounted() {
    this.updateSidebarDimensions()
    window.addEventListener('scroll', this.updateSidebarDimensions)
    window.addEventListener('resize', this.updateSidebarDimensions)
  },
  watch: {
    locale(newLocale, oldLocale) {
      if (newLocale !== oldLocale) {
        // 语言切换时清除当前标签，避免标签不匹配问题
        if (this.currentTag) {
          this.clearTag()
        } else {
          this.loadPostData()
        }
      }
    }
  },
  methods: {

    
    async loadPostData() {
      try {
        this.postData = await loadPosts() || [];
      } catch (error) {
        console.error('Failed to load post data:', error);
        this.postData = [];
      }
    },
    updateSidebarDimensions() {
      const header = document.querySelector('header')
      const footer = document.querySelector('footer')
      const content = this.$refs.sidebarContent
      const sidebarContainer = this.$refs.sidebarContainer

      if (!content || !sidebarContainer) return

      const headerHeight = header?.offsetHeight || 0
      const footerHeight = footer?.offsetHeight || 0
      const viewportHeight = window.innerHeight
      const scrollTop = window.scrollY
      const documentHeight = document.documentElement.scrollHeight

      const remainingPageHeight = Math.max(
        0,
        documentHeight - scrollTop - headerHeight - footerHeight - 40
      )
      const availableHeight = Math.min(
        viewportHeight - headerHeight - 40,
        remainingPageHeight
      )

      content.style.maxHeight = `${availableHeight}px`
      content.style.overflowY = content.scrollHeight > availableHeight ? 'auto' : 'visible'
    },
    clearTag() {
      const q = { ...this.$route.query }
      delete q.tag
      q.page = '1'
      this.$router.push({ path: this.$route.path, query: q }).catch(() => {})
      this.$nextTick(() => {
        window.scrollTo({ top: 0, behavior: 'smooth' })
        // 重新加载数据以确保状态一致
        this.loadPostData()
      })
    }
  },
  beforeUnmount() {
    window.removeEventListener('scroll', this.updateSidebarDimensions)
    window.removeEventListener('resize', this.updateSidebarDimensions)
  }
}
</script>

<style scoped>
.sticky-sidebar {
  position: sticky;
  top: 30px;
  box-sizing: border-box;
  width: 100%;
  -webkit-overflow-scrolling: touch;
  transition: max-height 0.2s ease;
}

@media (max-width: 991px) {
  .sticky-sidebar {
    position: static;
    top: auto;
    bottom: auto !important;
    max-height: none !important;
    overflow-y: visible !important;
  }
}
.current-tag-chip {
  font-size: 1rem; /* 与卡片标签更一致 */
  font-weight: 500;
  color: var(--app-chip-text);
  background: var(--app-chip-bg); /* 白色背底 */
  padding: 0.25rem 0.5rem; /* 适度内边距 */
  border: 1px solid var(--app-chip-border); /* 细边框 */
  box-shadow: var(--app-card-shadow); /* 轻微卡片感 */
  border-radius: 12px; /* 圆角矩形 */
  align-items: center; /* 垂直居中 */
  gap: 0.4rem;
}

.chip-close {
  font-size: 1.2rem; /* 与芯片字体更协调 */
  line-height: 1;
  background: transparent;
  border: none;
  color: var(--app-chip-close-text); /* 与全局 text-secondary 接近 */
  padding: 0;
  margin-left: 2px; /* 与文字间距更自然 */
  cursor: pointer;
  opacity: 1; /* 常驻显示 */
  transition: color 0.2s ease;
}

.chip-close:hover {
  color: var(--app-primary); /* 与全局主色一致的悬停色 */
}
</style>
