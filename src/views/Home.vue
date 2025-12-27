<!-- Home.vue -->
<template>
  <div class="container view-container home-view">
    <div class="row py-4 px-0">
      <div class="col-12 col-lg-9 order-1 order-lg-2 typography-body mb-4 mb-lg-0" ref="mainContent">
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

/*
  HomeView
  - 家页面
*/
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
      if (window.innerWidth < 992) return

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
      this.$router.push({ path: this.$route.path, query: q }).catch(() => { })
      this.$nextTick(() => {
        window.scrollTo({ top: 0, behavior: 'smooth' })
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

.current-tag-chip {
  font-size: 1rem;
  font-weight: 500;
  color: var(--app-chip-text);
  background: var(--app-chip-bg);
  padding: 0.25rem 0.5rem;
  border: 1px solid var(--app-chip-border);
  box-shadow: var(--app-card-shadow);
  border-radius: 12px;
  align-items: center;
  gap: 0.4rem;
}

.chip-close {
  font-size: 1.2rem;
  line-height: 1;
  background: transparent;
  border: none;
  color: var(--app-chip-close-text);
  padding: 0;
  margin-left: 2px;
  cursor: pointer;
  opacity: 1;
  transition: color 0.2s ease;
}

.chip-close:hover {
  color: var(--app-primary);
}
</style>
